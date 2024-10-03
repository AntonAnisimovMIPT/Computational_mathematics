use nalgebra::{DMatrix, DVector};

pub fn bicgstab_solve(matrix: &DMatrix<f64>, f: &DVector<f64>, tolerance: f64, max_iterations: usize) -> (DVector<f64>, Vec<f64>) {
    let n = matrix.nrows();
    let mut x = DVector::zeros(n);  
    let mut r = f - matrix * &x;     
    let r_hat = r.clone();       
    let mut residuals = vec![];      
    let mut rho_old = 1.0;
    let mut alpha = 1.0;
    let mut omega = 1.0;
    let mut v = DVector::zeros(n);
    let mut p = DVector::zeros(n);

    for iteration in 0..max_iterations {
        let rho_new = r_hat.dot(&r);
        if rho_new.abs() < 1e-50 {
            break;
        }

        if iteration == 0 {
            p = r.clone();
        } else {
            let beta = (rho_new / rho_old) * (alpha / omega);
            p = &r + beta * (&p - omega * &v);
        }

        v = matrix * &p;
        alpha = rho_new / r_hat.dot(&v);
        let s = &r - alpha * &v;

        let residual_norm = s.norm();
        residuals.push(residual_norm);
        if residual_norm < tolerance {
            println!("BiCGSTAB method converged in {} iterations", iteration + 1);
            x += alpha * &p;
            return (x, residuals);
        }

        let t = matrix * &s;
        omega = t.dot(&s) / t.dot(&t);
        x += alpha * &p + omega * &s;
        r = &s - omega * &t;

        let residual_norm = r.norm();
        residuals.push(residual_norm);
        if residual_norm < tolerance {
            println!("BiCGSTAB method converged in {} iterations", iteration + 1);
            return (x, residuals);
        }

        rho_old = rho_new;
    }

    println!("Warning!!! BiCGSTAB method did not converge in {} iterations.", max_iterations);
    (x, residuals)
}
