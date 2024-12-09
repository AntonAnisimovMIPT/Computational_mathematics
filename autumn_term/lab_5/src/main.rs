fn trapezoidal_method_h(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() - 1;
    let mut integral = 0.0;

    for i in 0..n {
        let h = x[i + 1] - x[i];
        integral += h * (y[i] + y[i + 1]) / 2.0;
    }

    integral
}

fn trapezoidal_method_2h(x: &[f64], y: &[f64]) -> f64 {
    let mut integral = 0.0;

    for i in (0..x.len() - 1).step_by(2) {
        let h = x[i + 2] - x[i];
        integral += h * (y[i] + y[i + 2]) / 2.0;
    }

    integral
}

fn simpson_method(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() - 1;
    if n % 2 != 0 {
        panic!("Simpson's rule requires an even number of intervals.");
    }

    let h = (x[n] - x[0]) / n as f64;
    let mut integral = y[0] + y[n]; 

    for i in 1..n {
        if i % 2 == 0 {
            integral += 2.0 * y[i]; 
        } else {
            integral += 4.0 * y[i]; 
        }
    }

    integral * h / 3.0
}

fn runge_refinement(i_h: f64, i_2h: f64, p: u32) -> f64 {
    let r: f64 = 2.0; 
    (r.powi(p as i32) * i_h - i_2h) / (r.powi(p as i32) - 1.0)
}

fn main() {

    let x = vec![0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    let y = vec![
        0.0, 0.021470, 0.293050, 0.494105, 0.541341, 0.516855, 0.468617, 0.416531, 0.367879,
    ];

    let integral_trapezoidal = trapezoidal_method_h(&x, &y);
    let integral_trapezoidal_fine = trapezoidal_method_2h(&x, &y);

    let refined_result = runge_refinement(integral_trapezoidal, integral_trapezoidal_fine, 2);

    let integral_simpson = simpson_method(&x, &y);

    println!("Интеграл методом трапеций: {:.6}", integral_trapezoidal);
    println!(
        "Уточнённый интеграл методом Рунге: {:.6}",
        refined_result
    );
    println!("Интеграл методом Симпсона: {:.6}", integral_simpson);
}