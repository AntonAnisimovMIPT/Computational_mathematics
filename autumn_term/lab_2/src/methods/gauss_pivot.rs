use nalgebra::{DMatrix, DVector};

pub fn gauss_pivot_solve(matrix: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    let n = matrix.nrows();
    let mut a = matrix.clone();
    let mut b = b.clone();

    for i in 0..n {
        let max_row = (i..n).max_by(|&x, &y| a[(x, i)].abs().partial_cmp(&a[(y, i)].abs()).unwrap()).unwrap();

        if max_row != i {
            a.swap_rows(i, max_row);
            b.swap_rows(i, max_row);
        }

        for j in i + 1..n {
            let factor = a[(j, i)] / a[(i, i)];
            for k in i..n {
                a[(j, k)] -= factor * a[(i, k)];
            }
            b[j] -= factor * b[i];
        }
    }

    let mut x = DVector::zeros(n);
    for i in (0..n).rev() {
        let mut sum = 0.0;
        for j in i + 1..n {
            sum += a[(i, j)] * x[j];
        }
        x[i] = (b[i] - sum) / a[(i, i)];
    }

    x
}

