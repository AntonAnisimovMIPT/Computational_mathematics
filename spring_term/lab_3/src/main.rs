use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufWriter, Write};

const H: f64 = 0.005;
const N: usize = (1.0 / H) as usize + 1;

fn p_squared(x: f64) -> f64 {
    10.0 + (2.0 * PI * x).sin()
}

fn f_func(x: f64) -> f64 {
    (2.0 * PI * x).cos()
}

fn solve_periodic() -> Vec<f64> {
    let x: Vec<f64> = (0..N).map(|i| i as f64 * H).collect();

    let a: Vec<f64> = vec![1.0; N];
    let b: Vec<f64> = x
        .iter()
        .map(|&xi| 2.0 + p_squared(xi) * H.powi(2))
        .collect();
    let c: Vec<f64> = vec![1.0; N];
    let d: Vec<f64> = x.iter().map(|&xi| f_func(xi) * H.powi(2)).collect();

    let mut alpha = vec![0.0; N];
    let mut beta = vec![0.0; N];
    let mut gamma = vec![0.0; N];

    alpha[1] = c[0] / b[0];
    beta[1] = -d[0] / b[0];
    gamma[1] = a[0] / b[0];

    for n in 1..N - 1 {
        let den = b[n] - a[n] * alpha[n];
        alpha[n + 1] = c[n] / den;
        beta[n + 1] = (-d[n] + a[n] * beta[n]) / den;
        gamma[n + 1] = (a[n] * gamma[n]) / den;
    }

    let den_N = -b[N - 1] + a[N - 1] * (alpha[N - 1] + gamma[N - 1]);
    let mu_N = -c[N - 1] / den_N;
    let nu_N = (d[N - 1] - a[N - 1] * beta[N - 1]) / den_N;

    let mut mu = vec![0.0; N];
    let mut nu = vec![0.0; N];
    mu[N - 1] = mu_N;
    nu[N - 1] = nu_N;

    for n in (1..N).rev() {
        mu[n - 1] = alpha[n] * mu[n] + gamma[n] * mu_N;
        nu[n - 1] = beta[n] + alpha[n] * nu[n] + gamma[n] * nu_N;
    }

    let y0 = nu[0] / (1.0 - mu[0]);
    let mut y = vec![0.0; N];
    y[0] = y0;
    y[N - 1] = mu[N - 1] * y[0] + nu[N - 1];

    for n in (1..N).rev() {
        y[n - 1] = alpha[n] * y[n] + beta[n] + gamma[n] * y[N - 1];
    }

    y
}

fn save_results(x: &[f64], y: &[f64]) {
    let file = File::create("data.txt").expect("Unable to create file");
    let mut writer = BufWriter::new(file);

    for (xi, yi) in x.iter().zip(y.iter()) {
        writeln!(writer, "{:.5} {:.10}", xi, yi).expect("Write error");
    }
}

fn main() {
    let y = solve_periodic();
    let x: Vec<f64> = (0..N).map(|i| i as f64 * H).collect();
    save_results(&x, &y);
}
