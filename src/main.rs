use aero::Panel;
use fa::prelude::*;
use fa::FaerMat;
use sample_points::sample_points;

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
const N: usize = 4;
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

fn get_vel_at(
    at: (f64, f64),
    sol: &faer::Mat<f64>,
    panels: &Vec<Panel>,
    freestream: (f64, f64),
) -> (f64, f64) {
    let mut vel = freestream;
    //let mut vel = (0.0, 0.0);
    let vortex_strength = sol[(panels.len(), 0)];
    for (idx, panel) in panels.iter().enumerate() {
        let source_strength = sol[(idx, 0)];
        let source_vel = panel.source_vel_at(at);
        let vortex_vel = panel.vortex_vel_at(at);

        vel.0 += source_vel.0 * source_strength + vortex_vel.0 * vortex_strength;
        vel.1 += source_vel.1 * source_strength + vortex_vel.1 * vortex_strength;
        //vel.0 += source_vel.0;
        //vel.1 += source_vel.1;
        //vel.0 += vortex_vel.0;
        //vel.1 += vortex_vel.1;
    }

    vel
}

fn vel_field_panel(
    sol: &faer::Mat<f64>,
    panels: &Vec<Panel>,
    freestream: (f64, f64),
) -> Vec<((f64, f64), (f64, f64))> {
    let points = panels.iter().map(|panel| panel.midpoint());

    points
        .into_iter()
        .filter_map(|point| {
            let mut vel = get_vel_at(point, sol, panels, freestream);
            vel.0 *= 0.1;
            vel.1 *= 0.1;
            if vel.0 * vel.0 + vel.1 * vel.1 > 5.0 {
                None
            } else {
                Some((point, vel))
            }
        })
        .collect()
}

fn vel_field_solved(
    sol: &faer::Mat<f64>,
    panels: &Vec<Panel>,
    freestream: (f64, f64),
) -> Vec<((f64, f64), (f64, f64))> {
    let points = sample_points((-0.2, -0.1), (1.2, 0.25));

    points
        .into_iter()
        .filter_map(|point| {
            let mut vel = get_vel_at(point, sol, panels, freestream);
            vel.0 *= 0.1;
            vel.1 *= 0.1;
            if vel.0 * vel.0 + vel.1 * vel.1 > 5.0 {
                None
            } else {
                Some((point, vel))
            }
        })
        .collect()
}

fn main() {
    let points = get_points();
    let mut min = (0.0, 0.0);
    let mut max = (0.0, 0.0);
    for point in &points {
        if point.0 < min.0 {
            min.0 = point.0
        }
        if point.1 < min.1 {
            min.1 = point.1
        }
        if point.0 > max.0 {
            max.0 = point.0
        }
        if point.1 > max.1 {
            max.1 = point.1
        }
    }
    let gmin = (min.0 - 0.1, min.1 - 0.1);
    let gmax = (max.0 + 0.1, max.1 + 0.1);
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
    if DO_PLOTS {
        /*let nrms: Vec<((f64, f64), (f64, f64))> = panels
            .iter()
            .map(|panel| (panel.midpoint(), panel.normal()))
            .collect();
        plots::vector_field::plot_vector_field("normals", nrms, (gmin, gmax), 0.1);*/
    }

    let vstrength = 0.2;
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
            if cause_idx == effect_idx {
                approx::assert_relative_eq!(source.0, 0.0, epsilon = 0.0000001);
                approx::assert_relative_eq!(source.1, 0.0, epsilon = 0.0000001);
                approx::assert_relative_eq!(vortex.0, 0.0, epsilon = 0.0000001);
                approx::assert_relative_eq!(vortex.1, 0.0, epsilon = 0.0000001);
            }

            // Source is scaled by sigma_cause (one for each panel)
            mat[(effect_idx, cause_idx)] = source_dot;
            // Vortex is scaled by "global" vortex intensity
            // (cause panels.len() represents the vortex)
            mat[(effect_idx, panels.len())] += vortex_dot;
        }

        // Right hand side
        rhs[(effect_idx, 0)] = -(freestream.0 * nrm.0 + freestream.1 * nrm.1);
    }

    // Kutta condition on trailing edge, which implies that flow is parallel to it
    // We may approximate it by considering the resulting velocity on the two
    // trailing edge panels.
    // TODO: Make sure both trailing edge panels are same length?
    /*let top = &panels[0];
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

        // TODO: Demonstrate that this is equivalent
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

            println!("{:?} {:?}", source_cond, vortex_cond);

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
    */
    mat[(panels.len(), panels.len())] = 1.0;
    rhs[(panels.len(), 0)] = vstrength;

    println!("{:?}", mat);
    println!("{:?}", rhs);

    let plu = mat.full_piv_lu();
    let x = plu.solve(&rhs);

    if DO_PLOTS {
        let vel_field = vel_field_solved(&x, &panels, freestream);
        plots::vector_field::plot_vector_field(
            "velfield",
            vel_field,
            ((-0.2, -0.1), (1.2, 0.25)),
            0.05,
        );

        let NUM_STREAMLINES = 150;
        let stream_start: Vec<(f64, f64)> = (0..NUM_STREAMLINES)
            .into_iter()
            .map(|i| (i as f64) / (NUM_STREAMLINES as f64))
            .map(|yprog| yprog * (gmax.1 - gmin.1) + gmin.1)
            .map(|y| (gmin.0, y))
            .collect();

        plots::stream_plot::plot_stream_plot(
            "stream",
            stream_start,
            |point| get_vel_at(point, &x, &panels, freestream),
            (gmin, gmax),
            0.01,
            1.5,
            Some(&points),
        );
    }

    println!("{:?}", x);

    for p in &panels {
        let cvel = get_vel_at(p.midpoint(), &x, &panels, freestream);
        let nrm0 = cvel.0 * p.normal().0;
        let nrm1 = cvel.1 * p.normal().1;
        assert_relative_eq!(nrm0, -nrm1, max_relative = 0.001);
    }
}
