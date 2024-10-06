use aero::Panel;

extern crate faer as fa;
extern crate std;

#[allow(unused_imports)]
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
const DO_TESTS: bool = false;

fn get_points() -> Vec<(f64, f64)> {
    let mut upper_resampled = geom::resample(&points::POINTS_UPPER, N_POINTS_UPPER);
    upper_resampled.reverse();
    upper_resampled.pop();
    let lower_resampled = geom::resample(&points::POINTS_LOWER, N_POINTS_LOWER);

    approx::assert_relative_eq!(upper_resampled[0].0, lower_resampled.last().unwrap().0);
    approx::assert_relative_eq!(upper_resampled[0].1, lower_resampled.last().unwrap().1);

    // combine both, erasing shared element

    [upper_resampled, lower_resampled].concat()
}

fn get_panels(points: &Vec<(f64, f64)>) -> Vec<Panel> {
    let mut out = Vec::<Panel>::new();

    for point_pair in points.windows(2) {
        let npanel = Panel {
            start: point_pair[0],
            end: point_pair[1],
        };
        out.push(npanel);
    }

    return out;
}

fn single_panel_cases() {
    visual_tests::single_panel::single_panel_horizontal();
    visual_tests::single_panel::single_panel_backwards();
    visual_tests::single_panel::single_panel_tilted();
    visual_tests::single_panel::single_panel_tilted_backwards();
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
        if DO_TESTS {
            single_panel_cases();
            compound_cases();
        }
        plots::geom::plot(&points);
    }

    // Build the panels
    let panels = get_panels(&points);

    // Each row of the matrix represents an equation, each column an unknown
    // thus a (effect, cause) maps to (row, column)
    let mut mat = fa::Mat::<f64>::zeros(panels.len() + 1, panels.len() + 1);
    for (effect_idx, effect) in panels.iter().enumerate() {
        // Each panel has at its center the equation:
        // (sum of induced velocities + freestream) * normal = 0
        let nrm = effect.normal();
        let midpoint = effect.midpoint();
        for (cause_idx, cause) in panels.iter().enumerate() {
            let source = cause.source_vel_at(midpoint);
            let vortex = cause.vortex_vel_at(midpoint);

            let source_dot = source.0 * nrm.0 + source.1 * nrm.1;
            let vortex_dot = vortex.0 * nrm.0 + vortex.1 * nrm.1;

            assert!(!f64::is_nan(source_dot));
            assert!(!f64::is_nan(vortex_dot));

            // Source is scaled by sigma_cause (one for each panel)
            mat[(effect_idx, cause_idx)] = source_dot;
            // Vortex is scaled by "global" vortex intensity
            mat[(panels.len(), cause_idx)] = vortex_dot;
        }
    }

    // Kutta condition on trailing edge
}
