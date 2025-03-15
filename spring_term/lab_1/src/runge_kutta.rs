pub fn runge_kutta1<F>(
    f: F,
    x0: f64,
    z0: f64,
    t0: f64,
    t_end: f64,
    h: f64,
    _e: f64,
) -> Vec<(f64, f64, f64)>
where
    F: Fn(f64, f64, f64) -> (f64, f64),
{
    let mut t = t0;
    let mut x = x0;
    let mut z = z0;
    let mut result = Vec::new();

    while t <= t_end {
        result.push((t, x, z));

        let (fx, fz) = f(x, z, t);
        x += h * fx;
        z += h * fz;

        t += h;
    }

    result
}

pub fn runge_kutta2<F>(
    f: F,
    x0: f64,
    z0: f64,
    t0: f64,
    t_end: f64,
    h: f64,
    _e: f64,
) -> Vec<(f64, f64, f64)>
where
    F: Fn(f64, f64, f64) -> (f64, f64),
{
    let mut t = t0;
    let mut x = x0;
    let mut z = z0;
    let mut result = Vec::new();

    while t <= t_end {
        result.push((t, x, z));

        let (k1x, k1z) = f(x, z, t);
        let (k2x, k2z) = f(x + 0.5 * h * k1x, z + 0.5 * h * k1z, t + 0.5 * h);

        x += h * k2x;
        z += h * k2z;

        t += h;
    }

    result
}

pub fn runge_kutta3<F>(
    f: F,
    x0: f64,
    z0: f64,
    t0: f64,
    t_end: f64,
    h: f64,
    _e: f64,
) -> Vec<(f64, f64, f64)>
where
    F: Fn(f64, f64, f64) -> (f64, f64),
{
    let mut t = t0;
    let mut x = x0;
    let mut z = z0;
    let mut result = Vec::new();

    while t <= t_end {
        result.push((t, x, z));

        let (k1x, k1z) = f(x, z, t);
        let (k2x, k2z) = f(x + 0.5 * h * k1x, z + 0.5 * h * k1z, t + 0.5 * h);
        let (k3x, k3z) = f(x - h * k1x + 2.0 * h * k2x, z - h * k1z + 2.0 * h * k2z, t + h);

        x += (h / 6.0) * (k1x + 4.0 * k2x + k3x);
        z += (h / 6.0) * (k1z + 4.0 * k2z + k3z);

        t += h;
    }

    result
}

pub fn runge_kutta4<F>(f: F, x0: f64, z0: f64, t0: f64, t_end: f64, h: f64, e: f64) -> Vec<(f64, f64, f64)>
where
    F: Fn(f64, f64, f64) -> (f64, f64),
{
    let mut t = t0;
    let mut x = x0;
    let mut z = z0;
    let mut result = Vec::new();

    while t <= t_end {
        let (k1x, k1z) = f(x, z, t);
        let (k2x, k2z) = f(x + 0.5 * h * k1x, z + 0.5 * h * k1z, t + 0.5 * h);
        let (k3x, k3z) = f(x + 0.5 * h * k2x, z + 0.5 * h * k2z, t + 0.5 * h);
        let (k4x, k4z) = f(x + h * k3x, z + h * k3z, t + h);

        x += (h / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
        z += (h / 6.0) * (k1z + 2.0 * k2z + 2.0 * k3z + k4z);

        t += h;
        result.push((t, x, z));
    }

    result
}
