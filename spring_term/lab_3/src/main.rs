use std::fs::File;
use std::io::{Write, BufWriter};
use std::f64::consts::PI;

const H: f64 = 0.005;
const N: usize = (1.0 / H) as usize;

fn p2(x: f64) -> f64 {
    10.0 + (2.0 * PI * x).sin()
}

fn f(x: f64) -> f64 {
    (2.0 * PI * x).cos()
}

fn solve_periodic() -> Vec<f64> {
    let a = vec![-1.0; N];
    let mut b = vec![0.0; N];
    let c = vec![-1.0; N];
    let mut d = vec![0.0; N];
    
    for i in 0..N {
        let x = i as f64 * H;
        b[i] = 2.0 + H * H * p2(x);
        d[i] = H * H * f(x);
    }
    
    let mut alpha = vec![0.0; N];
    let mut beta = vec![0.0; N];
    
    alpha[0] = c[0] / b[0];
    beta[0] = d[0] / b[0];
    
    for i in 1..N-1 {
        let denom = b[i] - a[i] * alpha[i - 1];
        alpha[i] = c[i] / denom;
        beta[i] = (d[i] - a[i] * beta[i - 1]) / denom;
    }
    
    let denom = b[N-1] - a[N-1] * alpha[N-2] - c[N-1] * alpha[0];
    let beta_n = (d[N-1] - a[N-1] * beta[N-2] - c[N-1] * beta[0]) / denom;
    
    let mut y = vec![0.0; N];
    y[N-1] = beta_n;
    for i in (0..N-1).rev() {
        y[i] = beta[i] - alpha[i] * y[i + 1];
    }
    
    y
}

fn save_to_file(data: &Vec<f64>) {
    let file = File::create("data.txt").expect("Unable to create file");
    let mut writer = BufWriter::new(file);
    
    for (i, &val) in data.iter().enumerate() {
        writeln!(writer, "{:.5} {:.5}", i as f64 * H, val).expect("Write error");
    }
}

fn main() {
    let solution = solve_periodic();
    save_to_file(&solution);
}