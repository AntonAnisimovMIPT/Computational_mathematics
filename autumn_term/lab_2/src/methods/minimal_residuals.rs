use nalgebra::{DMatrix, DVector};

pub fn minimal_residuals_solve(matrix: &DMatrix<f64>, f: &DVector<f64>, tolerance: f64, max_iterations: usize) -> (DVector<f64>, Vec<f64>) {
    let n = matrix.nrows();
    let mut x = DVector::zeros(n);  
    let mut r = f - matrix * &x;     
    let mut residuals = vec![];     

    for iteration in 0..max_iterations {
        let z = matrix * &r;  
        let alpha = r.dot(&r) / r.dot(&z);  
        x += alpha * &r;  
        r = f - matrix * &x;  

        let residual_norm = r.norm();  
        residuals.push(residual_norm);

        if residual_norm < tolerance {
            println!("Minimal residual method converged in {} iterations", iteration + 1);
            return (x, residuals);
        }
    }

    println!("Warning!!! Minimal residual method did not converge in {} iterations.", max_iterations);
    (x, residuals)
}