use nalgebra::{DMatrix, DVector};

pub fn gradient_descent_solve(matrix: &DMatrix<f64>, f: &DVector<f64>, learning_rate: f64, tolerance: f64, max_iterations: usize) -> (DVector<f64>, Vec<f64>) {
    let n = matrix.nrows();
    let mut x = DVector::zeros(n); 
    let mut r = f - matrix * &x; 
    let mut residuals = Vec::new(); 

    for iteration in 0..max_iterations {
        
        let gradient = matrix.transpose() * &r;

        x += learning_rate * &gradient;

        r = f - matrix * &x;
        residuals.push(r.norm()); 

        if r.norm() < tolerance {
            println!("Gradient descent converged in {} iterations", iteration + 1);
            return (x, residuals);
        }
    }

    println!("Warning!!! Gradient descent did not converge in {} iterations.", max_iterations);
    (x, residuals)
}
