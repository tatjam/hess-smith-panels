use fa::prelude::*;
use fa::FaerMat;

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

pub struct SplineSegment {
    t0: f64,
    t1: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
}

pub fn sample_spline(spline: &Vec<SplineSegment>, x: f64) -> f64 {
    for seg in spline.iter() {
        if seg.t0 <= x && seg.t1 >= x {
            return x.powi(3) * seg.a + x.powi(2) * seg.b + x.powi(1) * seg.c + seg.d;
        }
    }

    assert!(false, "Segment outside spline");
    return 0.0;
}

// Generates pts.len() - 1 spline segments, already solved
// If rounding != 0, instead of imposing y'' = 0 at the last point, it imposes
// y' = rounding
pub fn find_spline(pts: &[(f64, f64)], rounding: f64) -> Vec<SplineSegment> {
    assert!(pts.len() > 1);

    let mut out = Vec::new();

    // Build linear system
    // (Each spline segment has 4 equations and 4 unknowns)
    let mut mat = fa::Mat::zeros(4 * (pts.len() - 1), 4 * (pts.len() - 1));
    let mut rhs = fa::Mat::zeros(4 * (pts.len() - 1), 1);

    if (rounding == 0.0) {
        // Extreme segment, left second derivative null
        // (1 equation)
        mat[(0, 0)] = 6.0 * pts[0].0.powi(1);
        mat[(0, 1)] = 2.0;
        rhs[(0, 0)] = 0.0;
    } else {
        // Extreme segment, first derivative set by rounding
        // (1 equation)
        mat[(0, 0)] = 3.0 * pts[0].0.powi(2);
        mat[(0, 1)] = 2.0 * pts[0].0.powi(1);
        mat[(0, 2)] = 1.0;
        rhs[(0, 0)] = rounding;
    }
    // Central segments
    for i in 0..(pts.len() - 1) {
        let eq_offset = if i == 0 {
            1
        } else {
            // Zeroth segment only had 3 equations
            i * 4 - 1
        };
        let unknown_offset = i * 4;

        // First equation: aX^3 + bX^2 + cX + d = Y
        // Left continuity equation
        mat[((eq_offset + 0), (unknown_offset + 0))] = pts[i].0.powi(3);
        mat[((eq_offset + 0), (unknown_offset + 1))] = pts[i].0.powi(2);
        mat[((eq_offset + 0), (unknown_offset + 2))] = pts[i].0.powi(1);
        mat[((eq_offset + 0), (unknown_offset + 3))] = 1.0;
        rhs[(eq_offset, 0)] = pts[i].1;

        // Second equation a(X_next)^3 + b(X_next)^2 + c(X_next) + d = Y_next
        // Right continuity equation
        // Not present on right-most point
        let doff = if i != pts.len() - 1 {
            mat[((eq_offset + 1), (unknown_offset + 0))] = pts[i + 1].0.powi(3);
            mat[((eq_offset + 1), (unknown_offset + 1))] = pts[i + 1].0.powi(2);
            mat[((eq_offset + 1), (unknown_offset + 2))] = pts[i + 1].0.powi(1);
            mat[((eq_offset + 1), (unknown_offset + 3))] = 1.0;
            rhs[(eq_offset + 1, 0)] = pts[i + 1].1;
            1
        } else {
            0
        };

        if i != 0 {
            let prev_eq_unknown_offset = (i - 1) * 4;
            // Third equation (first deriv continuity at the left)
            // 3a_prev X^2 + 2b_prev X + c_prev = 3aX^2 + 2bX + c
            // (Grouping terms on left side, as they are all non constants)
            mat[((eq_offset + 1 + doff), (prev_eq_unknown_offset + 0))] = 3.0 * pts[i].0.powi(2);
            mat[((eq_offset + 1 + doff), (prev_eq_unknown_offset + 1))] = 2.0 * pts[i].0.powi(1);
            mat[((eq_offset + 1 + doff), (prev_eq_unknown_offset + 2))] = 1.0;
            mat[((eq_offset + 1 + doff), (unknown_offset + 0))] = -3.0 * pts[i].0.powi(2);
            mat[((eq_offset + 1 + doff), (unknown_offset + 1))] = -2.0 * pts[i].0.powi(1);
            mat[((eq_offset + 1 + doff), (unknown_offset + 2))] = -1.0;
            rhs[(eq_offset + 1 + doff, 0)] = 0.0;

            // Fourth equation (second deriv continuity at the left)
            // 6a_prev X + 2b_prev = 6aX + 2b
            // (Grouping terms on left side, as they are all non constants)
            mat[((eq_offset + 2 + doff), (prev_eq_unknown_offset + 0))] = 6.0 * pts[i].0.powi(1);
            mat[((eq_offset + 2 + doff), (prev_eq_unknown_offset + 1))] = 2.0;
            mat[((eq_offset + 2 + doff), (unknown_offset + 0))] = -6.0 * pts[i].0.powi(1);
            mat[((eq_offset + 2 + doff), (unknown_offset + 1))] = -2.0;
            rhs[(eq_offset + 2 + doff, 0)] = 0.0;
        }
    }

    // Extreme segment, right second derivative null
    // 1 equation
    let last_eq_offset = (pts.len() - 1) * 4 - 1;
    let last_eq_unknown_offset = (pts.len() - 2) * 4;

    mat[(last_eq_offset, last_eq_unknown_offset + 0)] = 6.0 * pts[pts.len() - 1].0.powi(1);
    mat[(last_eq_offset, last_eq_unknown_offset + 1)] = 2.0;
    rhs[(last_eq_offset, 0)] = 0.0;

    let sol = mat.full_piv_lu().solve(rhs);

    for i in 0..(pts.len() - 1) {
        let seg = SplineSegment {
            t0: pts[i].0,
            t1: pts[i + 1].0,
            a: sol[(i * 4 + 0, 0)],
            b: sol[(i * 4 + 1, 0)],
            c: sol[(i * 4 + 2, 0)],
            d: sol[(i * 4 + 3, 0)],
        };
        out.push(seg);
    }

    return out;
}

/// Cosine sampling is used, to accumulate points near the leading and trialing edges
/// desired_n is the number of points to be placed
/// Always keeps first and last points the same
pub fn resample(
    pts: &[(f64, f64)],
    desired_n: usize,
    use_spline: bool,
    rounding: f64,
) -> Vec<(f64, f64)> {
    assert!(is_monotonically_increasing_in_x(pts));

    let (min, max) = find_x_extremes(pts);
    let mut out: Vec<(f64, f64)> = Vec::new();

    if use_spline {
        let spline = find_spline(pts, rounding);
        for n in 0..desired_n {
            let prog = (n as f64) / (desired_n as f64 - 1.0);
            let theta = prog * std::f64::consts::PI; // Theta goes from 0->PI
            let cos_sample = 0.5 - 0.5 * theta.cos();
            let x = cos_sample * (max - min) + min;
            let y = sample_spline(&spline, x);

            out.push((x, y));
        }
    } else {
        for n in 0..desired_n {
            let prog = (n as f64) / (desired_n as f64 - 1.0);
            let theta = prog * std::f64::consts::PI; // Theta goes from 0->PI
            let cos_sample = 0.5 - 0.5 * theta.cos();
            let x = cos_sample * (max - min) + min;
            let ((x0, y0), (x1, y1)) =
                find_containing_segment(pts, x).expect("Profile data malformed");

            // Linear interpolation within the segment
            let in_segment_progress = (x - x0) / (x1 - x0);
            let y = (y1 - y0) * in_segment_progress + y0;

            out.push((x, y));
        }
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
        let resampled1 = resample(&[(0.0, 0.0), (1.0, 1.0)], N, false, 0.0);
        let resampled2 = resample(&[(0.0, 0.0), (0.5, 0.5), (1.0, 1.0)], N, false, 0.0);
        let resampled3 = resample(
            &[(0.0, 0.0), (0.1, 0.1), (0.65, 0.65), (1.0, 1.0)],
            N,
            false,
            0.0,
        );

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

    #[test]
    fn test_spline() {
        let pts = [(0.0, 0.0), (0.5, 0.5), (1.0, 1.0)];
        let spline = find_spline(&pts, 0.0);
        assert_relative_eq!(sample_spline(&spline, 0.5), 0.5);
    }
}
