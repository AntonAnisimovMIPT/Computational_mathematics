pub struct CubicSpline {
    a: Vec<f64>,  
    b: Vec<f64>,  
    c: Vec<f64>,  
    d: Vec<f64>,  
    x: Vec<f64>,  
}

impl CubicSpline {
    pub fn new(x: &Vec<f64>, y: &Vec<f64>) -> Self {
        let n = x.len();
        let a = y.clone();
        let mut b = vec![0.0; n - 1];
        let mut d = vec![0.0; n - 1];
        let mut h = vec![0.0; n - 1];
        let mut alpha = vec![0.0; n - 1];
        
        for i in 0..n - 1 {
            h[i] = x[i + 1] - x[i];
            alpha[i] = if i > 0 {
                (3.0 / h[i]) * (a[i + 1] - a[i]) - (3.0 / h[i - 1]) * (a[i] - a[i - 1])
            } else {
                0.0
            };
        }

        let mut c = vec![0.0; n];
        let mut l = vec![1.0; n];
        let mut mu = vec![0.0; n];
        let mut z = vec![0.0; n];
        
        for i in 1..n - 1 {
            l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        for j in (0..n - 1).rev() {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }

        CubicSpline { a, b, c, d, x: x.to_vec() }
    }

    pub fn evaluate(&self, x: f64) -> f64 {
        let n = self.x.len();
        let mut i = n - 2;

        for j in 0..n - 1 {
            if x >= self.x[j] && x <= self.x[j + 1] {
                i = j;
                break;
            }
        }

        let dx = x - self.x[i];
        self.a[i] + self.b[i] * dx + self.c[i] * dx.powi(2) + self.d[i] * dx.powi(3)
    }
}
