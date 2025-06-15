mod lightcurve;
mod config;
mod util;

use std::env;
use std::process;
use rayon::prelude::*;
use std::fs::OpenOptions;
use std::io::Write;
use lightcurve::Lightcurve;
use config::Config;
use rayon::prelude::*;
use std::f64::consts::FRAC_PI_2;
// use std::io::{BufRead};
// use rand::Rng;

fn main() {
    let args: Vec<String> = env::args().collect();

    let config = Config::build(&args).unwrap_or_else(|err| {
        println!("Problem parsing arguments: {err}");
        process::exit(1);
    });

    let lightcurves = Lightcurve::load_lightcurves(
        config.input,
        config.redshifts,
        &config.mjd_name,
        &config.mag_name,
        &config.magerr_name
    );

    let (tau, sf, num) = calc_esf(
        &lightcurves,
        &config.tau_lo,
        &config.tau_hi
    );

    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(config.output).unwrap();

    let line = tau.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
    writeln!(file, "{}", line).unwrap();

    let line = sf.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
    writeln!(file, "{}", line).unwrap();

    let line = num.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
    writeln!(file, "{}", line).unwrap();

    file.flush().unwrap();

    if !config.bootstrap {
        return;
    }

    for _ in 0..1000 {
        let sample = Lightcurve::resample(&lightcurves);
        let (_, sf, _) = calc_esf(
            &sample,
            &config.tau_lo,
            &config.tau_hi
        );

        let line = sf.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        writeln!(file, "{}", line).unwrap();
    }
}


fn parallel_sum_all_mag_diffs(
    lightcurves: &Vec<Lightcurve>,
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
    lightcurves: &Vec<Lightcurve>,
    tau_lo: &[f64],
    tau_hi: &[f64]
) -> (Vec<f64>, Vec<f64>, Vec<usize>) {
    let mut res_tau = Vec::<f64>::new();
    let mut res_sf = Vec::<f64>::new();
    let mut res_num = Vec::<usize>::new();

    for (lo, hi) in tau_lo.iter().zip(tau_hi.iter()) {
        let (sum_tau, sum_diff, sum_diff_err_square, n) =
            parallel_sum_all_mag_diffs(lightcurves, *lo, *hi);

        let tau = sum_tau / n as f64;
        let diff = sum_diff / n as f64;
        let diff_err_square = sum_diff_err_square / n as f64;

        let sf = (FRAC_PI_2 * diff * diff - diff_err_square).sqrt();
        res_tau.push(tau);
        res_sf.push(sf);
        res_num.push(n);
    }

    (res_tau, res_sf, res_num)
}
