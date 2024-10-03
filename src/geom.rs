/// Returns (min, max) pair
fn find_x_extremes(pts: &[(f64, f64)]) -> (f64, f64) {
    pts.iter()
        .fold((999999999.0, -999999999.0), |(min, max), (x, _)| {
            let nmin = if *x < min { *x } else { min };
            let nmax = if *x > max { *x } else { max };
            (nmin, nmax)
        })
}

fn is_monotonically_increasing_in_x(pts: &[(f64, f64)]) -> bool {
    let mut prev_x: f64 = -999999999999.0;
    for (x, _) in pts {
        if *x < prev_x {
            return false;
        }
        prev_x = *x;
    }
    return true;
}

// Returns (start_point, end_point) of containing segment
pub fn find_containing_segment(pts: &[(f64, f64)], x: f64) -> Option<((f64, f64), (f64, f64))> {
    for sl in pts.windows(2) {
        let (x0, y0) = sl[0];
        let (x1, y1) = sl[1];
        if x0 <= x && x1 >= x {
            return Some(((x0, y0), (x1, y1)));
        }
    }
    return None;
}

/// Cosine sampling is used, to accumulate points near the leading and trialing edges
/// desired_n is the number of points to be placed
/// Always keeps first and last points the same
pub fn resample(pts: &[(f64, f64)], desired_n: usize) -> Vec<(f64, f64)> {
    assert!(is_monotonically_increasing_in_x(pts));

    let (min, max) = find_x_extremes(pts);
    let mut out: Vec<(f64, f64)> = Vec::new();

    for n in 0..desired_n {
        let prog = (n as f64) / (desired_n as f64 - 1.0);
        let theta = prog * std::f64::consts::PI; // Theta goes from 0->PI
        let cos_sample = 0.5 - 0.5 * theta.cos();
        let x = cos_sample * (max - min) + min;
        let ((x0, y0), (x1, y1)) = find_containing_segment(pts, x).expect("Profile data malformed");

        // Linear interpolation within the segment
        let in_segment_progress = (x - x0) / (x1 - x0);
        let y = (y1 - y0) * in_segment_progress + y0;

        out.push((x, y));
    }

    return out;
}

#[cfg(test)]
mod test {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn test_find_x_extremes() {
        let (min, max) = find_x_extremes(&[
            (132.0, 3.0),
            (-2.0, 1.0),
            (-43.0, 3.0),
            (137.0, 3.0),
            (137.0, -543.0),
        ]);
        assert_eq!(min, -43.0);
        assert_eq!(max, 137.0);
    }

    #[test]
    fn test_resample() {
        use std::f64::*;

        const N: usize = 10;
        let resampled1 = resample(&[(0.0, 0.0), (1.0, 1.0)], N);
        let resampled2 = resample(&[(0.0, 0.0), (0.5, 0.5), (1.0, 1.0)], N);
        let resampled3 = resample(&[(0.0, 0.0), (0.1, 0.1), (0.65, 0.65), (1.0, 1.0)], N);

        assert_eq!(resampled1.len(), N);
        assert_eq!(resampled2.len(), N);
        assert_eq!(resampled3.len(), N);

        for t in 0..N {
            let theta = (t as f64) / (N as f64 - 1.0);
            let v = 0.5 * (1.0 - f64::cos(consts::PI * theta));
            assert_relative_eq!(resampled1[t].0, v);
            assert_relative_eq!(resampled1[t].1, v);
            assert_relative_eq!(resampled2[t].0, v);
            assert_relative_eq!(resampled2[t].1, v);
            assert_relative_eq!(resampled3[t].0, v);
            assert_relative_eq!(resampled3[t].1, v);
        }
    }
}
