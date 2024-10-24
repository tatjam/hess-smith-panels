const NUM_POINTS: isize = 20;

pub fn sample_points(start: (f64, f64), end: (f64, f64)) -> Vec<(f64, f64)> {
    assert!(start.0 < end.0);
    assert!(start.1 < end.1);

    let mut out = Vec::<(f64, f64)>::new();

    for xi in -NUM_POINTS..NUM_POINTS {
        for yi in -NUM_POINTS..NUM_POINTS {
            let x = (xi as f64) / (NUM_POINTS as f64) * (end.0 - start.0) + start.0;
            let y = (yi as f64) / (NUM_POINTS as f64) * (end.1 - start.1) + start.1;
            out.push((x, y));
        }
    }
    return out;
}
