use std::env;
use std::process;
use sfrs::Config;

fn main() {
    let args: Vec<String> = env::args().collect();

    let config = Config::build(&args).unwrap_or_else(|err| {
        println!("Problem parsing arguments: {err}");
        process::exit(1);
    });

    let lightcurves = sfrs::Lightcurve::load_lightcurves(
        config.input,
        config.redshifts,
        &config.mjd_name,
        &config.mag_name,
        &config.magerr_name
    );

    let (tau, sf) = sfrs::calc_esf(
        &lightcurves,
        &config.tau_lo,
        &config.tau_hi
    );

    println!("{:?}", tau);
    println!("{:?}", sf);
}
