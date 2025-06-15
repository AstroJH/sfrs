use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use super::util;

pub struct Config {
    pub input: Vec<String>,
    pub output: String,
    pub redshifts: Vec<f64>,
    pub tau_lo: Vec<f64>,
    pub tau_hi: Vec<f64>,
    pub mjd_name: String,
    pub mag_name: String,
    pub magerr_name: String,
    pub bootstrap: bool
}

impl Config {
    pub fn build(args: &[String]) -> Result<Config, &'static str> {
        println!("{:?}", args);
        if args.len() < 5 {
            return Err("not enough arguments");
        }

        let input_file = args[1].clone();
        let output = args[2].clone();
        let redshift_file = args[3].clone();
        let taubin_file = args[4].clone();

        let rtd = env::var("SFRS_RTD").unwrap_or_else(|_| ".".to_string());
        let bootstrap = env::var("SFRS_BOOTSTRAP").is_ok();
        let mjd_name = env::var("SFRS_MJD").unwrap_or_else(|_| "mjd".to_string());
        let mag_name = env::var("SFRS_MAG").unwrap_or_else(|_| "mag".to_string());
        let magerr_name = env::var("SFRS_MAGERR").unwrap_or_else(|_| "magerr".to_string());

        let (tau_lo, tau_hi) = Config::load_taubins(&taubin_file);
        let redshifts = Config::load_nums(&redshift_file);
        let input = Config::load_inputs(&input_file, &rtd);


        Ok(Config {
            input,
            output,
            redshifts,
            tau_lo,
            tau_hi,
            mjd_name,
            mag_name,
            magerr_name,
            bootstrap
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

    fn build_taubins(method: &str, start: f64, end: f64, num: usize) -> (Vec<f64>, Vec<f64>) {
        let result = match method {
            "linspace" => {
                Ok(util::linspace(start, end, num, true))
            }

            "logspace" => {
                Ok(util::logspace(start, end, num, 10., true))
            }

            "geomspace" => {
                Ok(util::geomspace(start, end, num, true))
            }

            _ => Err("Invalid method.")
        }.unwrap();

        (
            Vec::<f64>::from(&result[..result.len() - 1]),
            Vec::<f64>::from(&result[1..])
        )
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
    use super::*;
    #[test]
    fn test_build_taubins() {
        let (tau_lo, tau_hi) =
            Config::build_taubins("logspace", 1., 5., 5);

        assert_eq!(tau_lo, vec![10.,  100.,  1000.,  10000.]);
        assert_eq!(tau_hi, vec![100., 1000., 10000., 100000.]);
    }
}