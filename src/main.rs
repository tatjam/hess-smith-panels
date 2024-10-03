extern crate faer as fa;
extern crate std;

#[macro_use]
extern crate approx;

mod geom;
mod plots;
mod points;

const N: usize = 70;
const N_POINTS_UPPER: usize = N / 2 + 1;
const N_POINTS_LOWER: usize = N / 2 + 1;

fn get_points() -> Vec<(f64, f64)> {
    let mut upper_resampled = geom::resample(&points::POINTS_UPPER, N_POINTS_UPPER);
    upper_resampled.reverse();
    let mut lower_resampled = geom::resample(&points::POINTS_LOWER, N_POINTS_LOWER);

    approx::assert_relative_eq!(upper_resampled[0].0, lower_resampled.last().unwrap().0);
    approx::assert_relative_eq!(upper_resampled[0].1, lower_resampled.last().unwrap().1);

    // combine both, erasing shared element
    lower_resampled.pop();

    [upper_resampled, lower_resampled].concat()
}

fn main() {
    let points = get_points();
    assert_eq!(points.len(), N + 1);

    plots::geom::plot(&points);
}
