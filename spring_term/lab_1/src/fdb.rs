use crate::runge_kutta;
pub fn fdb<F>(f: F, x0: f64, z0: f64, t0: f64, t_end: f64, h: f64, e: f64) -> Vec<(f64, f64, f64)>
where
    F: Fn(f64, f64, f64) -> (f64, f64),
{
    let mut t = t0;
    let mut x = x0;
    let mut z = z0;
    let mut result = Vec::new();

    let initial = runge_kutta(&f, x0, z0, t0, t0 + 2.0 * h, h, e);
    result.extend_from_slice(&initial);

    while t <= t_end {
        let (fx, fz) = f(x, z, t);
        let (fx_prev, fz_prev) = f(x - h, z - h, t - h);

        let dx = h * (fx + fx_prev) / 2.0;
        let dz = h * (fz + fz_prev) / 2.0;

        x += dx;
        z += dz;

        t += h;
        result.push((t, x, z));
    }

    result
}
