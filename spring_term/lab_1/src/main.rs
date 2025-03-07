use std::fs::File;
use std::io::Write;

mod runge_kutta;
mod adams;
mod fdb;
use crate::{adams::adams, fdb::fdb, runge_kutta::runge_kutta};

fn main() {
    let e = 0.1;
    let h = 0.01;
    let t0 = 0.0;
    let t_end = 100.0;
    let x0 = 2.0;
    let z0 = 0.0;

    let f = |x: f64, z: f64, _t: f64| (z, e * (1.0 - x * x) * z - x);

    let result = runge_kutta(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/runge_kutta.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("runge_kutta computated");

    let result = adams(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/adams.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("adams computated");

    let result = fdb(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/fdb.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("fdb computated");
}