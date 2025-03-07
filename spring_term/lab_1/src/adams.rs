use crate::runge_kutta;
pub fn adams<F>(f: F, x0: f64, z0: f64, t0: f64, t_end: f64, h: f64, e: f64) -> Vec<(f64, f64, f64)>
where
    F: Fn(f64, f64, f64) -> (f64, f64),
{
    let mut t = t0;
    let mut x = x0;
    let mut z = z0;
    let mut result = Vec::new();

    let initial = runge_kutta(&f, x0, z0, t0, t0 + 2.0 * h, h, e);
    result.extend_from_slice(&initial);

    let mut x_prev = initial[1].1;
    let mut z_prev = initial[1].2;

    while t <= t_end {
        let (fx, fz) = f(x, z, t);
        let (fx_prev, fz_prev) = f(x_prev, z_prev, t - h);
        let (fx_prev2, fz_prev2) = f(x_prev - h, z_prev - h, t - 2.0 * h);

        let dx = h * (3.0 / 2.0 * fx - 1.0 / 2.0 * fx_prev);
        let dz = h * (3.0 / 2.0 * fz - 1.0 / 2.0 * fz_prev);

        x += dx;
        z += dz;
        t += h;

        result.push((t, x, z));

        x_prev = x;
        z_prev = z;
    }

    result
}
