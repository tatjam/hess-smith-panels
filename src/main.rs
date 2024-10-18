use aero::Panel;
use fa::prelude::*;
use fa::FaerMat;
use sample_points::sample_points;

extern crate plotters as pl;
use plotters::prelude::*;

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
const N: usize = 100;
const N_POINTS_UPPER: usize = N / 2 + 1;
const N_POINTS_LOWER: usize = N / 2 + 1;
const DO_PLOTS: bool = true;
const DO_TESTS: bool = false;

fn get_points() -> Vec<(f64, f64)> {
    let mut upper_resampled = geom::resample(&points::POINTS_UPPER, N_POINTS_UPPER, true);
    upper_resampled.reverse();
    upper_resampled.pop();
    let lower_resampled = geom::resample(&points::POINTS_LOWER, N_POINTS_LOWER, true);

    //approx::assert_relative_eq!(upper_resampled[0].0, lower_resampled.last().unwrap().0);
    //approx::assert_relative_eq!(upper_resampled[0].1, lower_resampled.last().unwrap().1);

    // combine both, erasing shared element at leading edge

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

        let svel = panel.source_vel_at(at);
        let tvel = panel.vortex_vel_at(at);

        vel.0 += source_strength * svel.0 + vortex_strength * tvel.0;
        vel.1 += source_strength * svel.1 + vortex_strength * tvel.1;
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
            if vel.0 * vel.0 + vel.1 * vel.1 > 5.0 {
                None
            } else {
                Some((point, vel))
            }
        })
        .collect()
}

fn find_plot_extremes(points: &Vec<(f64, f64)>) -> ((f64, f64), (f64, f64)) {
    let mut min = (0.0, 0.0);
    let mut max = (0.0, 0.0);
    for point in points {
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
    return (gmin, gmax);
}

fn main() {
    let points = get_points();
    assert_eq!(points.len(), N + 1);

    let (gmin, gmax) = find_plot_extremes(&points);

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
        let nrms: Vec<((f64, f64), (f64, f64))> = panels
            .iter()
            .map(|panel| (panel.midpoint(), panel.normal()))
            .collect();
        plots::vector_field::plot_vector_field("normals", nrms, (gmin, gmax), 0.1);
    }

    // Each row of the matrix represents an equation, each column an unknown
    // thus a (effect, cause) maps to (row, column)
    let u_infty: f64 = 1.0;
    let alpha: f64 = 0.0;
    let freestream = (u_infty * alpha.cos(), u_infty * alpha.sin());
    let mut mat = fa::Mat::<f64>::zeros(panels.len() + 1, panels.len() + 1);
    let mut rhs = fa::Mat::<f64>::zeros(panels.len() + 1, 1);
    for (effect_idx, effect) in panels.iter().enumerate() {
        let theta_i = effect.theta();
        let midpoint = effect.midpoint();

        for (cause_idx, cause) in panels.iter().enumerate() {
            let theta_j = cause.theta();
            let params = cause.params_at(midpoint);
            // Source is scaled by sigma_cause (one for each panel)
            mat[(effect_idx, cause_idx)] = 0.5
                * std::f64::consts::FRAC_1_PI
                * (params.0 * (theta_i - theta_j).sin() + params.1 * (theta_i - theta_j).cos());
            // Vortex is scaled by "global" vortex intensity
            // (cause panels.len() represents the vortex)
            mat[(effect_idx, panels.len())] += 0.5
                * std::f64::consts::FRAC_1_PI
                * (-params.1 * (theta_i - theta_j).sin() + params.0 * (theta_i - theta_j).cos())
                * params.2;
        }

        // Right hand side
        rhs[(effect_idx, 0)] = -u_infty * (alpha - theta_i).sin();
    }

    // Kutta condition
    // TODO: Check this or the other way around
    let top = panels.last().unwrap();
    let bottom = panels.first().unwrap();
    let theta_n = top.theta();
    let theta_1 = bottom.theta();
    for (cause_idx, cause) in panels.iter().enumerate() {
        let params_n = cause.params_at(top.midpoint());
        let params_1 = cause.params_at(bottom.midpoint());
        let theta_j = cause.theta();

        mat[(panels.len(), cause_idx)] = 0.5
            * std::f64::consts::FRAC_1_PI
            * (params_1.1 * (theta_1 - theta_j).sin() - params_1.0 * (theta_1 - theta_j).cos()
                + params_n.1 * (theta_n - theta_j).sin()
                - params_n.0 * (theta_n - theta_j).cos());

        mat[(panels.len(), panels.len())] += 0.5
            * std::f64::consts::FRAC_1_PI
            * (params_1.1 * (theta_1 - theta_j).cos() + params_1.0 * (theta_1 - theta_j).sin())
            * params_1.2;
        mat[(panels.len(), panels.len())] += 0.5
            * std::f64::consts::FRAC_1_PI
            * (params_n.1 * (theta_n - theta_j).cos() + params_n.0 * (theta_n - theta_j).sin())
            * params_n.2;
    }

    rhs[(panels.len(), 0)] = -u_infty * ((alpha - theta_1).cos() + (alpha - theta_n).cos());

    let plu = mat.full_piv_lu();
    let x = plu.solve(&rhs);

    if DO_PLOTS {
        let vel_field = vel_field_panel(&x, &panels, freestream);
        plots::vector_field::plot_vector_field(
            "velfield",
            vel_field,
            ((-0.2, -0.1), (1.2, 0.25)),
            0.05,
        );

        let NUM_STREAMLINES = 50;
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

    let mut cps: Vec<(f64, f64)> = Vec::new();
    for p in &panels {
        let params = p.params_at((0.0, 0.1));
        let cvel = get_vel_at(p.midpoint(), &x, &panels, freestream);
        let nrm0 = cvel.0 * p.normal().0;
        let nrm1 = cvel.1 * p.normal().1;
        assert_relative_eq!(nrm0, -nrm1, max_relative = 0.001);

        let tvelmag2 = cvel.0 * cvel.0 + cvel.1 * cvel.1;

        let cpi = 1.0 - tvelmag2 / u_infty.powi(2);
        cps.push((p.midpoint().0, -cpi));
    }

    let root = BitMapBackend::new("out/cps.png", (1024, 1024)).into_drawing_area();

    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .margin(15)
        .set_left_and_bottom_label_area_size(35)
        .build_cartesian_2d(-0.2..1.2, -2.0..1.5)
        .unwrap();

    chart.draw_series(LineSeries::new(cps, &RED)).unwrap();

    chart.configure_mesh().draw().unwrap();
}
