use std::fs::File;
use std::sync::{Arc};
use csv::Reader;
use rayon::prelude::*;
use rand::Rng;

#[derive(Debug, Clone)]
struct LightcurveData {
    mjd: Box<[f64]>,
    mag: Box<[f64]>,
    magerr: Box<[f64]>,
}

#[derive(Debug, Clone)]
pub struct Lightcurve {
    data: Arc<LightcurveData>, // Read-only
    pub redshift: f64
}

impl Lightcurve {
    pub fn new(
        mjd: Vec<f64>,
        mag: Vec<f64>,
        magerr: Vec<f64>,
        redshift: f64
    ) -> Self {
        assert_eq!(mjd.len(), mag.len(), "mjd and mag must have same length");
        assert_eq!(mjd.len(), magerr.len(), "mjd and magerr must have same length");

        Lightcurve {
            data: Arc::new(LightcurveData {
                mjd: mjd.into_boxed_slice(),
                mag: mag.into_boxed_slice(),
                magerr: magerr.into_boxed_slice()
            }),
            redshift,
        }
    }

    pub fn len(&self) -> usize {
        self.data.mjd.len()
    }

    pub fn load(
        filename: &str, redshift: f64,
        mjd_name: &str, mag_name: &str, magerr_name: &str
    ) -> Result<Self, &'static str> {
        let file = File::open(filename);
        if file.is_err() {
            return Err("Error opening file.");
        }
        let file = file.unwrap();

        let mut reader = Reader::from_reader(file);

        let headers = reader.headers().unwrap();
        let (mut i_mjd, mut i_mag, mut i_magerr) = (0usize, 0usize, 0usize);

        let mut flag = 0u8;
        for (i, header) in headers.iter().enumerate() {
            if header.eq(mjd_name) {
                i_mjd = i;
                flag += 1;
            } else if header.eq(mag_name) {
                i_mag = i;
                flag += 1;
            } else if header.eq(magerr_name) {
                i_magerr = i;
                flag += 1;
            }

            if flag >=3 { break; }
        }

        if flag < 3 {
            panic!("ERROR!");
        }

        const CAPACITY: usize = 32;
        let mut mjd_vec = Vec::<f64>::with_capacity(CAPACITY);
        let mut mag_vec = Vec::<f64>::with_capacity(CAPACITY);
        let mut magerr_vec = Vec::<f64>::with_capacity(CAPACITY);

        for record in reader.records().into_iter() {
            let record = record.unwrap();
            let mjd = record.get(i_mjd).unwrap().parse::<f64>().unwrap();
            let mag = record.get(i_mag).unwrap().parse::<f64>().unwrap();
            let magerr = record.get(i_magerr).unwrap().parse::<f64>().unwrap();

            mjd_vec.push(mjd);
            mag_vec.push(mag);
            magerr_vec.push(magerr);
        }

        Ok(Lightcurve::new(mjd_vec, mag_vec, magerr_vec, redshift))
    }


    pub fn load_lightcurves(input: Vec<String>, redshift: Vec<f64>,
                            mjd_name: &str, mag_name: &str, magerr_name: &str) -> Vec<Self> {
        input.par_iter().zip(redshift.par_iter())
            .filter_map(|(filename, z)| {
                match Lightcurve::load(filename, *z, mjd_name, mag_name, magerr_name) {
                    Ok(lightcurve) => Some(lightcurve),
                    Err(_) => None
                }
            })
            .collect()
    }

    pub fn find_mag_differences(
        &self,
        tau_lo: f64,
        tau_hi: f64
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let mjd = &self.data.mjd;
        let mag = &self.data.mag;
        let magerr = &self.data.magerr;
        let n = mjd.len();

        if n < 2 { return (Vec::new(), Vec::new(), Vec::new()) }

        // 计算相邻时间间隔
        let mut intervals = Vec::with_capacity(n-1);
        for i in 1..mjd.len() {
            intervals.push(
                (mjd[i] - mjd[i-1])/(1.0+self.redshift)
            );
        }

        let mut mag_diffs = Vec::<f64>::new();
        let mut mag_diff_errs = Vec::<f64>::new();
        let mut taus = Vec::<f64>::new();

        // lp := left pointer
        // rp := right pointer
        for lp in 0..(n-1) {
            let mut cumulative = intervals[lp];
            let mut rp = lp+1;

            while rp < n-1 && cumulative < tau_lo {
                cumulative += intervals[rp];
                rp += 1;
            }

            if cumulative < tau_lo {
                // 不再有时间间隔 >= tau_lo, stop!
                break;
            } else if cumulative > tau_hi {
                // 当前 right pointer 已不再有时间间隔位于 [tau_lo, tau_hi]
                // 下一轮尝试
                continue;
            }

            // 开始收集本轮满足要求的时间间隔点对
            while rp < n {
                let diff = (mag[lp] - mag[rp]).abs();
                let tau = mjd[rp] - mjd[lp];
                let diff_err = magerr[lp].powi(2) + magerr[rp].powi(2);

                mag_diffs.push(diff);
                taus.push(tau);
                mag_diff_errs.push(diff_err);

                if rp < n-1 {
                    cumulative += intervals[rp];
                    if cumulative > tau_hi {
                        break;
                    }
                }
                rp += 1;
            }
        }

        (taus, mag_diffs, mag_diff_errs)
    }

    pub fn resample(original: &[Self]) -> Vec<Self> {
        let n = original.len();
        let mut rng = rand::rng();

        (0..n)
            .map(|_| {
                let idx = rng.random_range(0..n);
                let lc = &original[idx];

                lc.clone()
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use crate::{Lightcurve};

    fn load_a_sample() -> Result<Lightcurve, &'static str> {
        let home = std::env::var("HOME").unwrap();
        let rtd = PathBuf::from(home)
            .join("Repository")
            .join("Paliya_2024_Seyfert1")
            .join("ztflc_repro");

        let filepath = rtd.join("10667-58163-0639.csv");

        let lcurve = Lightcurve::load(
            filepath.to_str().unwrap(),
            0.3,
            "mjd",
            "mag",
            "magerr"
        );

        lcurve
    }

    #[test]
    fn test_ref() {
        let lcurve = load_a_sample().unwrap();
        let lc_clone = lcurve.clone();

        println!("{:p}", lcurve.data.mag);
        println!("{:p}", lc_clone.data.mag);
    }

    #[test]
    fn it_works() {
        let lcurve = load_a_sample().unwrap();
        println!("{:#?}", lcurve)
    }
}