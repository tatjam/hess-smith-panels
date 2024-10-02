extern crate faer as fa;
extern crate std;
#[macro_use]
extern crate approx;

mod points;

/// Returns (min, max) pair
fn find_x_extremes(pts: &[(f64, f64)]) -> (f64, f64) {
    pts.iter()
        .fold((999999999.0, -999999999.0), |(min, max), (x, _)| {
            let nmin = if *x < min { *x } else { min };
            let nmax = if *x > max { *x } else { max };
            (nmin, nmax)
        })
}

/// Cosine sampling is used, to accumulate points near the leading and trialing edges
/// desired_n is the number of points to be placed
/// Always keeps first and last points the same
fn resample(pts: &[(f64, f64)], desired_n: usize) -> Vec<(f64, f64)> {
    let out: Vec<(f64, f64)> = Vec::new();

    for n in 0..desired_n {}

    return out;
}

fn main() {
    println!("Hello, world!");
}

#[cfg(test)]
mod test {
    use approx::relative_eq;

    use super::*;

    #[test]
    fn test_x_extremes() {
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
            relative_eq!(resampled1[t].0, v);
            relative_eq!(resampled1[t].1, v);
            relative_eq!(resampled2[t].0, v);
            relative_eq!(resampled2[t].1, v);
            relative_eq!(resampled3[t].0, v);
            relative_eq!(resampled3[t].1, v);
        }
    }
}
