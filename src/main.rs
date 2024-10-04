extern crate faer as fa;
extern crate std;

#[macro_use]
extern crate approx;

mod aero;
mod geom;
mod plots;
mod points;
mod sample_points;
mod visual_tests;

// Number of PANELS, point number is N + 1
const N: usize = 70;
const N_POINTS_UPPER: usize = N / 2 + 1;
const N_POINTS_LOWER: usize = N / 2 + 1;
const DO_PLOTS: bool = true;

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

fn single_panel_cases() {
    visual_tests::single_panel::single_panel_horizontal();
    visual_tests::single_panel::single_panel_tilted();
    visual_tests::single_panel::single_panel_backwards();
    visual_tests::single_panel::single_panel_vertical();
    visual_tests::single_panel::single_panel_vertical_backwards();
}

fn compound_cases() {
    visual_tests::compound_panel::compound_rhombus();
    visual_tests::compound_panel::compound_square();
}

fn main() {
    let points = get_points();
    assert_eq!(points.len(), N + 1);

    if DO_PLOTS {
        //single_panel_cases();
        compound_cases();
        plots::geom::plot(&points);
    }
}
