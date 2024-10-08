use aero::Panel;
use fa::prelude::*;
use fa::FaerMat;

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
    let freestream = (1.0, 0.0);
    let mut mat = fa::Mat::<f64>::zeros(panels.len() + 1, panels.len() + 1);
    let mut rhs = fa::Mat::<f64>::zeros(panels.len() + 1, 1);
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

        // Right hand side
        rhs[(effect_idx, 0)] = -(freestream.0 * nrm.0 + freestream.1 * nrm.1);
    }

    // Kutta condition on trailing edge, which implies that flow is parallel to it
    // We may approximate it by considering the resulting velocity on the two
    // trailing edge panels.
    // TODO: Make sure both trailing edge panels are same length?
    let top = &panels[0];
    let bottom = &panels[panels.len() - 1];
    let top_mp = top.midpoint();
    let bottom_mp = bottom.midpoint();
    let top_vec = (
        (top.end.0 - top.start.0) / top.len(),
        (top.end.1 - top.start.1) / top.len(),
    );
    let bottom_vec = (
        (bottom.start.0 - bottom.end.0) / bottom.len(),
        (bottom.start.1 - bottom.end.1) / bottom.len(),
    );
    for (cause_idx, cause) in panels.iter().enumerate() {
        let v_source_top = cause.source_vel_at(top_mp);
        let v_source_bottom = cause.source_vel_at(bottom_mp);
        let v_vortex_top = cause.vortex_vel_at(top_mp);
        let v_vortex_bottom = cause.vortex_vel_at(bottom_mp);

        let source_cond = v_source_top.0 * top_vec.0
            + v_source_top.1 * top_vec.1
            + v_source_bottom.0 * bottom_vec.0
            + v_source_bottom.1 * bottom_vec.1;

        let vortex_cond = v_vortex_top.0 * top_vec.0
            + v_vortex_top.1 * top_vec.1
            + v_vortex_bottom.0 * bottom_vec.0
            + v_vortex_bottom.1 * bottom_vec.1;

        // Source condition for each panel
        mat[(panels.len(), cause_idx)] = source_cond;
        // Vortex condition summed, cummulative effect of all vortices
        mat[(panels.len(), panels.len())] += vortex_cond;
    }
    // Freestream term, passed to right hand side, thus negative signs
    rhs[(panels.len(), 0)] = -top_vec.0 * freestream.0
        - top_vec.1 * freestream.1
        - bottom_vec.0 * freestream.0
        - bottom_vec.1 * freestream.1;

    let plu = mat.partial_piv_lu();
    let x = plu.solve(&rhs);

    println!("{:?}", x);
}
