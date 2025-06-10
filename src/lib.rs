use rayon::prelude::*;
use std::f64::consts::FRAC_PI_2;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use csv::{Reader};

#[derive(Debug)]
pub struct Lightcurve {
    pub mjd: Vec<f64>,
    pub mag: Vec<f64>,
    pub magerr: Vec<f64>,
    pub redshift: f64
}

pub struct Config {
    pub input: Vec<String>,
    pub redshifts: Vec<f64>,
    pub tau_lo: Vec<f64>,
    pub tau_hi: Vec<f64>,
    pub mjd_name: String,
    pub mag_name: String,
    pub magerr_name: String
}


impl Lightcurve {
    pub fn load(filename: &str, redshift: f64,
                mjd_name: &str, mag_name: &str, magerr_name: &str) -> Result<Self, &'static str> {
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

        let mut mjd_vec = Vec::<f64>::new();
        let mut mag_vec = Vec::<f64>::new();
        let mut magerr_vec = Vec::<f64>::new();

        for record in reader.records() {
            let record = record.unwrap();
            let mjd = record.get(i_mjd).unwrap().parse::<f64>().unwrap();
            let mag = record.get(i_mag).unwrap().parse::<f64>().unwrap();
            let magerr = record.get(i_magerr).unwrap().parse::<f64>().unwrap();

            mjd_vec.push(mjd);
            mag_vec.push(mag);
            magerr_vec.push(magerr);
        }


        Ok(Lightcurve {
            mjd: mjd_vec,
            mag: mag_vec,
            magerr: magerr_vec,
            redshift
        })
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

    fn find_mag_differences(
        &self,
        tau_lo: f64,
        tau_hi: f64
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let mjd = &self.mjd;
        let mag = &self.mag;
        let magerr = &self.magerr;
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
}

fn parallel_sum_all_mag_diffs(
    lightcurves: &[Lightcurve],
    tau_lo: f64,
    tau_hi: f64
) -> (f64, f64, f64, usize) {
    lightcurves.par_iter();
    let (sum_tau, sum_diff, sum_diff_err_square, n): (f64, f64, f64, usize) =
        lightcurves.par_iter().map(|lc| {

            let (tau, diff, diff_err_square) =
                lc.find_mag_differences(tau_lo, tau_hi);
            let sum_tau = tau.iter().sum::<f64>();
            let sum_diff = diff.iter().sum::<f64>();
            let sum_diff_err_square = diff_err_square.iter().sum::<f64>();
            let n = tau.len();

            (sum_tau, sum_diff, sum_diff_err_square, n)
        }).reduce(
            || (0.0, 0.0, 0.0, 0),
            |(a1, a2, a3, a4), (b1, b2, b3, b4)| (a1+b1, a2+b2, a3+b3, a4+b4)
        );

    (sum_tau, sum_diff, sum_diff_err_square, n)
}


pub fn calc_esf(
    lightcurves: &[Lightcurve],
    tau_lo: &[f64],
    tau_hi: &[f64]
) -> (Vec<f64>, Vec<f64>) {
    let mut res_tau = Vec::<f64>::new();
    let mut res_sf = Vec::<f64>::new();

    for (lo, hi) in tau_lo.iter().zip(tau_hi.iter()) {
        let (sum_tau, sum_diff, sum_diff_err_square, n) =
            parallel_sum_all_mag_diffs(lightcurves, *lo, *hi);

        let tau = sum_tau / n as f64;
        let diff = sum_diff / n as f64;
        let diff_err_square = sum_diff_err_square / n as f64;

        let sf = (FRAC_PI_2 * diff * diff - diff_err_square).sqrt();
        res_tau.push(tau);
        res_sf.push(sf);
    }

    (res_tau, res_sf)
}


impl Config {
    pub fn build(args: &[String]) -> Result<Config, &'static str> {
        if args.len() < 8 {
            return Err("not enough arguments");
        }

        let rtd = args[1].clone();
        let input_file = args[2].clone();
        let redshift_file = args[3].clone();
        let taubin_file = args[4].clone();

        let mjd_name = args[5].clone();
        let mag_name = args[6].clone();
        let magerr_name = args[7].clone();

        let (tau_lo, tau_hi) = Config::load_taubins(&taubin_file);
        let redshifts = Config::load_nums(&redshift_file);
        let input = Config::load_inputs(&input_file, &rtd);


        Ok(Config {
            input,
            redshifts,
            tau_lo,
            tau_hi,
            mjd_name,
            mag_name,
            magerr_name
        })
    }

    fn load_taubins(taubin_file: &str) -> (Vec<f64>, Vec<f64>) {
        let file = File::open(taubin_file).unwrap();
        let reader = BufReader::new(file);

        let mut tau_lo = vec![];
        let mut tau_hi = vec![];

        let mut flag = 0u8;
        for line in reader.lines() {
            let nums = line.unwrap().split(',')
                .filter_map(|x| x.parse::<f64>().ok())
                .collect::<Vec<f64>>();

            if flag == 0 { tau_lo = nums; }
            else if flag == 1 { tau_hi = nums; }
            else { break; }

            flag += 1;
        }

        (tau_lo, tau_hi)
    }

    fn load_nums(file: &str) -> Vec<f64> {
        let file = File::open(file).unwrap();
        let reader = BufReader::new(file);

        let numbers: Vec<f64> = reader
            .lines()
            .filter_map(|line| line.ok())           // 过滤读取错误
            .filter_map(|line| line.trim().parse().ok()) // 过滤解析错误
            .collect();

        numbers
    }

    fn load_inputs(file: &str, rtd: &str) -> Vec<String> {
        let file = File::open(file).unwrap();
        let reader = BufReader::new(file);

        reader.lines().map(|line| {
            PathBuf::from(rtd).join(line.unwrap()).to_str().unwrap().to_string()
        }).collect()
    }
}


#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use crate::Lightcurve;

    #[test]
    fn it_works() {
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

        println!("{:#?}", lcurve)
    }
}
