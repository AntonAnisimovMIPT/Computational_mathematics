use nalgebra::{DMatrix, DVector};

fn lu_decomposition(matrix: &DMatrix<f64>) -> (DMatrix<f64>, DMatrix<f64>) {
    let n = matrix.nrows();
    let mut l = DMatrix::identity(n, n); 
    let mut u = matrix.clone(); 

    for i in 0..n {
        for j in (i + 1)..n {
            
            let factor = u[(j, i)] / u[(i, i)];
            l[(j, i)] = factor; 

            for k in i..n {
                u[(j, k)] -= factor * u[(i, k)];
            }
        }
    }

    (l, u)
}

fn solve_lower_triangular(l: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    let n = l.nrows();
    let mut y = DVector::zeros(n);

    for i in 0..n {
        
        let sum: f64 = (0..i).map(|k| l[(i, k)] * y[k]).sum();
        
        y[i] = b[i] - sum;
    }

    y
}

fn solve_upper_triangular(u: &DMatrix<f64>, y: &DVector<f64>) -> DVector<f64> {
    let n = u.nrows();
    let mut x = DVector::zeros(n);

    for i in (0..n).rev() {
        let sum: f64 = (i + 1..n).map(|j| u[(i, j)] * x[j]).sum();
        x[i] = (y[i] - sum) / u[(i, i)];
    }

    x
}

pub fn lu_solve(matrix: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    let (l, u) = lu_decomposition(matrix);

    //println!("L:\n{}", l);
    //println!("U:\n{}", u);

    let y = solve_lower_triangular(&l, b);
    let x = solve_upper_triangular(&u, &y);

    x
}

