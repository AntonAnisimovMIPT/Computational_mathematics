use std::fs::File;
use std::io::Write;

mod runge_kutta;
mod adams;
mod bdf;
use crate::{
    adams::{adams2, adams3, adams4}, 
    bdf::{bdf2, bdf3, bdf4}, 
    runge_kutta::{runge_kutta1, runge_kutta2, runge_kutta3, runge_kutta4}};

fn main() {
    let e = 0.1;
    let h = 0.001;
    let t0 = 0.0;
    let t_end = 100.0;
    let x0 = 2.0;
    let z0 = 0.0;

    let f = |x: f64, z: f64, _t: f64| (z, e * (1.0 - x * x) * z - x);

    /// методы Рунге-Кутты
    let result = runge_kutta1(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/runge_kutta1.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("runge_kutta1 computated");

    let result = runge_kutta2(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/runge_kutta2.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("runge_kutta2 computated");

    let result = runge_kutta3(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/runge_kutta3.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("runge_kutta3 computated");

    let result = runge_kutta4(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/runge_kutta4.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("runge_kutta4 computated");

    /// методы Адамса
    let result = adams2(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/adams2.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("adams2 computated");

    let result = adams3(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/adams3.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("adams3 computated");

    let result = adams4(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/adams4.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("adams4 computated");

    /// формула дифференцирования назад
    let result = bdf2(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/bdf2.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("bdf2 computated");

    let result = bdf3(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/bdf3.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("bdf3 computated");

    let result = bdf4(&f, x0, z0, t0, t_end, h, e);
    let mut file = File::create("./results/data/bdf4.txt").unwrap();
    for (t, x, z) in result {
        writeln!(file, "{} {} {}", t, x, z).unwrap();
    }
    println!("bdf4 computated");
}