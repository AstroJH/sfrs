use std::f64;

pub fn linspace(start: f64, end: f64, num: usize, endpoint: bool) -> Vec<f64> {
    if num == 0 {
        return Vec::new();
    } else if num == 1 {
        return vec![start];
    }

    let step = if endpoint {
        (end - start) / ((num - 1) as f64)
    } else {
        (end - start) / (num as f64)
    };

    (0..num)
        .map(|i| start + step * (i as f64))
        .collect()
}


pub fn logspace(start: f64, end: f64, num: usize, base: f64, endpoint: bool) -> Vec<f64> {
    let exponents = linspace(start, end, num, endpoint);

    exponents.iter().map(|&exp| base.powf(exp)).collect()
}


pub fn geomspace(start: f64, end: f64, num: usize, endpoint: bool) -> Vec<f64> {
    assert!(start > 0.0, "start must be positive");
    assert!(end > 0.0, "end must be positive");

    let log_start = start.log10();
    let log_end = end.log10();

    let log_values = linspace(log_start, log_end, num, endpoint);

    log_values.iter().map(|&log_val| 10.0_f64.powf(log_val)).collect()
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_logspace() {
        assert_eq!(
            logspace(0., 5., 5, 10., false),
            vec![1.0, 10.0, 100.0, 1000.0, 10000.0]
        );
    }
}