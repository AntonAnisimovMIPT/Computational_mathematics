use std::f64::consts::E;

fn trapezoidal_method(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() - 1;
    let mut integral = 0.0;

    for i in 0..n {
        let h = x[i + 1] - x[i];
        integral += h * (y[i] + y[i + 1]) / 2.0;
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

fn runge_refinement(trapezoidal_coarse: f64, trapezoidal_fine: f64) -> f64 {
    let p = 2.0; 
    trapezoidal_fine + (trapezoidal_fine - trapezoidal_coarse) / (2_f64.powf(p) - 1.0)
}

fn interpolate(x: &[f64], y: &[f64], xi: f64) -> f64 {
    for i in 0..x.len() - 1 {
        if xi >= x[i] && xi <= x[i + 1] {
            let t = (xi - x[i]) / (x[i + 1] - x[i]);
            return y[i] * (1.0 - t) + y[i + 1] * t;
        }
    }
    0.0
}



fn main() {
    
    let x = vec![
        0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    let y = vec![
        0.0, 0.021470, 0.293050, 0.494105, 0.541341, 0.516855, 0.468617, 0.416531, 0.367879,
    ];

    let integral_trapezoidal = trapezoidal_method(&x, &y);

    // Создание более плотной сетки для уточнения
    let x_fine: Vec<f64> = (0..=16).map(|i| i as f64 / 16.0).collect();
    let y_fine: Vec<f64> = x_fine.iter().map(|xi| interpolate(&x, &y, *xi)).collect();
    let integral_trapezoidal_fine = trapezoidal_method(&x_fine, &y_fine);
    let refined_result = runge_refinement(integral_trapezoidal, integral_trapezoidal_fine);

    let integral_simpson = simpson_method(&x, &y);

   
    println!("Интеграл методом трапеций: {:.6}", integral_trapezoidal);
    println!(
        "Уточненный интеграл методом Рунге: {:.6}",
        refined_result
    );
    println!("Интеграл методом Симпсона: {:.6}", integral_simpson);
}