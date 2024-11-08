fn f1(x: f64, y: f64) -> f64 {
    x.sin() - y - 1.32
}

fn f2(x: f64, y: f64) -> f64 {
    y.cos() - x + 0.85
}

pub fn print_sistem_func() {
    println!("Выбранная система уравнений:
    г). sinx - y = 1.32
        cosy - x = -0.85");
}

fn dfx(x: f64, y: f64, h: f64) -> f64 {
    (f1(x + h, y) - f1(x - h, y)) / (2.0 * h)
}

fn dfy(x: f64, y: f64, h: f64) -> f64 {
    (f1(x, y + h) - f1(x, y - h)) / (2.0 * h)
}

fn dfx2(x: f64, y: f64, h: f64) -> f64 {
    (f2(x + h, y) - f2(x - h, y)) / (2.0 * h)
}

fn dfy2(x: f64, y: f64, h: f64) -> f64 {
    (f2(x, y + h) - f2(x, y - h)) / (2.0 * h)
}

pub fn mpi_system(x0: f64, y0: f64, tol: f64, max_iter: usize) -> Option<(f64, f64)> {
    let (mut x, mut y) = (x0, y0);
    for _ in 0..max_iter {
        let x_new = y.cos() + 0.85;
        let y_new = x.sin() - 1.32;
        if (x_new - x).abs() < tol && (y_new - y).abs() < tol {
            if f1(x_new, y_new).abs() < tol && f2(x_new, y_new).abs() < tol {
                return Some((x_new, y_new));
            } else {
                println!(
                    "Ошибка: найденное решение методом МПИ не удовлетворяет системе с заданной точностью."
                );
                return None;
            }
        }
        x = x_new;
        y = y_new;
    }
    None
}

pub fn newton_system(x0: f64, y0: f64, tol: f64, max_iter: usize, h: f64) -> Option<(f64, f64)> {
    let (mut x, mut y) = (x0, y0);
    for _ in 0..max_iter {
        let fx = f1(x, y);
        let fy = f2(x, y);
        let j11 = dfx(x, y, h);
        let j12 = dfy(x, y, h);
        let j21 = dfx2(x, y, h);
        let j22 = dfy2(x, y, h);
        
        let det = j11 * j22 - j12 * j21;
        if det.abs() < tol { return None; } 
        
        let delta_x = (-fx * j22 + fy * j12) / det;
        let delta_y = (fx * j21 - fy * j11) / det;
        
        x += delta_x;
        y += delta_y;
        
        if delta_x.abs() < tol && delta_y.abs() < tol {
            if f1(x, y).abs() < tol && f2(x, y).abs() < tol {
                return Some((x, y));
            } else {
                println!(
                    "Ошибка: найденное решение методом Ньютона не удовлетворяет системе с заданной точностью."
                );
                return None;
            }
        }
    }
    None
}
