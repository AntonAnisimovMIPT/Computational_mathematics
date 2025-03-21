use std::f64::consts::LN_10;

fn f(x: f64) -> f64 {
    2.0 * x.ln() / LN_10 - x / 2.0 + 1.0
}

pub fn print_func() {
    println!("Выбранное уравнение: 
    з). 2lgx - x/2 + 1 = 0");
}

fn df(x: f64, h: f64) -> f64 {
    (f(x + h) - f(x - h)) / (2.0 * h)
}

pub fn mpi_1_with_errors(x0: f64, tolerance: f64, max_iter: usize) -> Option<(f64, Vec<f64>)> {
    let mut x = x0;
    let mut errors = Vec::new();
    for _ in 0..max_iter {
        let x_new = 10_f64.powf(x / 4.0 - 0.5);
        let error = (x_new - x).abs();
        errors.push(error);

        if error < tolerance {
            if f(x_new).abs() < tolerance {
                return Some((x_new, errors));
            } else {
                println!(
                    "Ошибка: корень, найденный методом МПИ (начальное приближение x0 = {}), не удовлетворяет точности по f(x): f(x) = {:.6}",
                    x0,
                    f(x_new)
                );
                return None;
            }
        }
        x = x_new;
    }
    None
}

pub fn mpi_2_with_errors(x0: f64, tolerance: f64, max_iter: usize) -> Option<(f64, Vec<f64>)> {
    let mut x = x0;
    let mut errors = Vec::new();
    for _ in 0..max_iter {
        let x_new = 4.0 * x.ln() / LN_10 + 2.0;
        let error = (x_new - x).abs();
        errors.push(error);

        if error < tolerance {
            if f(x_new).abs() < tolerance {
                return Some((x_new, errors));
            } else {
                println!(
                    "Ошибка: корень, найденный методом МПИ (начальное приближение x0 = {}), не удовлетворяет точности по f(x): f(x) = {:.6}",
                    x0,
                    f(x_new)
                );
                return None;
            }
        }
        x = x_new;
    }
    None
}

pub fn newton_with_errors(x0: f64, tolerance: f64, max_iter: usize, h: f64) -> Option<(f64, Vec<f64>)> {
    let mut x = x0;
    let mut errors = Vec::new();
    for _ in 0..max_iter {
        let fx = f(x);
        let dfx = df(x, h);
        if dfx.abs() < tolerance { return None; }
        let x_new = x - fx / dfx;
        let error = (x_new - x).abs();
        errors.push(error);

        if error < tolerance {
            if f(x_new).abs() < tolerance {
                return Some((x_new, errors));
            } else {
                println!(
                    "Ошибка: корень, найденный методом Ньютона (начальное приближение x0 = {}), не удовлетворяет точности по f(x): f(x) = {:.6}",
                    x0,
                    f(x_new)
                );
                return None;
            }
        }
        x = x_new;
    }
    None
}
