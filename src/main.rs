use aero::Panel;
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

const DO_TESTS: bool = false;
const USE_ROUNDING: bool = false;

const PLOT_MIN: (f64, f64) = (-0.1, -0.1);
const PLOT_MAX: (f64, f64) = (1.4, 0.2);

fn get_points(n: usize) -> Vec<(f64, f64)> {
    // Number of PANELS, point number is N + 1
    let n_points_upper: usize = n / 2 + 1;
    let n_points_lower: usize = n / 2 + 1;

    let mut upper_resampled = geom::resample(
        &points::POINTS_UPPER,
        n_points_upper,
        true,
        if USE_ROUNDING { -3.5 } else { 0.0 },
    );
    upper_resampled.reverse();
    upper_resampled.pop();
    let lower_resampled = geom::resample(
        &points::POINTS_LOWER,
        n_points_lower,
        true,
        if USE_ROUNDING { 3.5 } else { 0.0 },
    );

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

fn many_operating_points() {
    let points = get_points(200);
    assert_eq!(points.len(), 200 + 1);

    plots::geom::plot(&points);

    // Build the panels
    let panels = get_panels(&points);
    let nrms: Vec<((f64, f64), (f64, f64))> = panels
        .iter()
        .map(|panel| (panel.midpoint(), panel.normal()))
        .collect();
    plots::vector_field::plot_vector_field("normals", nrms, (PLOT_MIN, PLOT_MAX), 0.1);

    let u_infty: f64 = 1.0;

    // A bunch of operating points
    for alpha_deg in (-15.0..15.0).step(2.5).values() {
        let alpha = alpha_deg * std::f64::consts::PI / 180.0;
        operating_point(&panels, &points, u_infty, alpha);
    }
}

fn operating_point(panels: &[Panel], points: &[(f64, f64)], u_infty: f64, alpha: f64) {
    let freestream = (u_infty * alpha.cos(), u_infty * alpha.sin());

    let x = aero::hess_smith(panels, u_infty, alpha);

    let vel_field = aero::vel_field_panel(&x, panels, u_infty, alpha);
    let fname = format!("velfield-{}", alpha);
    plots::vector_field::plot_vector_field(&fname, vel_field, ((-0.2, -0.1), (1.2, 0.25)), 0.05);

    let NUM_STREAMLINES = 50;
    // NOTE: Streamlines flow backwards
    let stream_start: Vec<(f64, f64)> = (0..NUM_STREAMLINES)
        .into_iter()
        .map(|i| (i as f64) / (NUM_STREAMLINES as f64))
        .map(|yprog| PLOT_MIN.1 + yprog * (PLOT_MAX.1 - PLOT_MIN.1))
        .map(|y| (PLOT_MAX.0, y))
        .collect();

    let fname = format!("stream-{}", alpha);
    plots::stream_plot::plot_stream_plot(
        &fname,
        stream_start,
        |point| aero::get_vel_at(point, &x, &panels, u_infty, alpha),
        (PLOT_MIN, PLOT_MAX),
        0.01,
        2.0,
        Some(&points),
    );

    let mut cps: Vec<(f64, f64)> = Vec::new();
    for p in panels {
        let params = p.params_at((0.0, 0.1));
        let cvel = aero::get_vel_at(p.midpoint(), &x, &panels, u_infty, alpha);
        let nrm0 = cvel.0 * p.normal().0;
        let nrm1 = cvel.1 * p.normal().1;
        assert_relative_eq!(nrm0, -nrm1, max_relative = 0.001);

        let tvelmag2 = cvel.0 * cvel.0 + cvel.1 * cvel.1;

        let mut cpi = 1.0 - tvelmag2 / u_infty.powi(2);
        // Clamp very extreme values
        if cpi.abs() > 4.0 {
            cpi /= cpi.abs();
            cpi *= 4.0;
        }
        cps.push((p.midpoint().0, cpi));
    }

    let fname = format!("out/{}-cps.png", alpha);
    let root = BitMapBackend::new(&fname, (1024, 1024)).into_drawing_area();

    root.fill(&WHITE).unwrap();

    let min_cp = cps.iter().min_by(|a, b| a.1.total_cmp(&b.1)).unwrap().1;
    let max_cp = cps.iter().max_by(|a, b| a.1.total_cmp(&b.1)).unwrap().1;

    let mut chart = ChartBuilder::on(&root)
        .margin(15)
        .set_left_and_bottom_label_area_size(35)
        .build_cartesian_2d(-0.2..1.2, max_cp..min_cp)
        .unwrap();

    chart.draw_series(LineSeries::new(cps, &RED)).unwrap();

    let txt = format!("Alpha = {}deg", alpha * std::f64::consts::FRAC_1_PI * 180.0);
    let txt_style = TextStyle::from(("sans-serif", 40).into_font()).color(&BLACK);
    root.draw_text(&txt, &txt_style, (600, 200)).unwrap();

    chart
        .configure_mesh()
        .y_desc("cp")
        .x_desc("x (m)")
        .draw()
        .unwrap();
}

fn num_points_sweep() {}

fn alpha_sweep() {
    let points = get_points(200);
    assert_eq!(points.len(), 200 + 1);

    // Build the panels
    let panels = get_panels(&points);
    let u_infty: f64 = 1.0;
    // alpha, cL, cDi, cM in that order
    let mut points: Vec<(f64, f64, f64, f64)> = Vec::new();
    // A bunch of operating points
    for alpha_deg in (-30.0..30.0).step(5.0).values() {
        let alpha = alpha_deg * std::f64::consts::PI / 180.0;
        let x = aero::hess_smith(&panels, u_infty, alpha);

        // These coefficients are NOT divided by chord, but it's one!
        let (cl, cdi, cm) = aero::find_coeffs(&panels, &x, u_infty, alpha, (0.25, 0.0));

        points.push((alpha_deg, cl, cdi, cm));

        println!("alpha = {} cl={}, cdi={}, cm={}", alpha, cl, cdi, cm);
    }

    {
        let root = BitMapBackend::new("out/cl.png", (512, 512)).into_drawing_area();

        root.fill(&WHITE).unwrap();

        let mut chart = ChartBuilder::on(&root)
            .margin(15)
            .set_left_and_bottom_label_area_size(35)
            .build_cartesian_2d(-30.0..30.0, -4.0..4.0)
            .unwrap();

        chart
            .draw_series(LineSeries::new(points.iter().map(|p| (p.0, p.1)), &RED))
            .unwrap()
            .label("cL (hess-smith)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &RED));

        chart
            .draw_series(LineSeries::new(
                points.iter().map(|p| {
                    (
                        p.0,
                        2.0 * std::f64::consts::PI * p.0 * std::f64::consts::PI / 180.0 + 0.203635,
                    )
                }),
                &BLUE,
            ))
            .unwrap()
            .label("cL (2pi ideal)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &BLUE));

        chart
            .configure_mesh()
            .x_desc("alpha (deg)")
            .y_desc("cL")
            .draw()
            .unwrap();

        chart
            .configure_series_labels()
            .background_style(&WHITE)
            .draw()
            .unwrap()
    }

    {
        let root = BitMapBackend::new("out/cdi.png", (512, 512)).into_drawing_area();

        root.fill(&WHITE).unwrap();

        let mut chart = ChartBuilder::on(&root)
            .margin(15)
            .set_left_and_bottom_label_area_size(35)
            .build_cartesian_2d(-30.0..30.0, -0.005..0.012)
            .unwrap();

        chart
            .draw_series(LineSeries::new(points.iter().map(|p| (p.0, p.2)), &RED))
            .unwrap()
            .label("cDi (hess-smith)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &RED));

        chart
            .configure_mesh()
            .x_desc("alpha (deg)")
            .y_desc("cDi")
            .draw()
            .unwrap();
        chart
            .configure_series_labels()
            .background_style(&WHITE)
            .draw()
            .unwrap()
    }

    {
        let root = BitMapBackend::new("out/cm.png", (512, 512)).into_drawing_area();

        root.fill(&WHITE).unwrap();

        let mut chart = ChartBuilder::on(&root)
            .margin(15)
            .set_left_and_bottom_label_area_size(35)
            .build_cartesian_2d(-30.0..30.0, -0.12..0.05)
            .unwrap();

        chart
            .draw_series(LineSeries::new(points.iter().map(|p| (p.0, p.3)), &RED))
            .unwrap()
            .label("cM (hess-smith)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &RED));

        chart
            .configure_mesh()
            .x_desc("alpha (deg)")
            .y_desc("cM")
            .draw()
            .unwrap();
        chart
            .configure_series_labels()
            .background_style(&WHITE)
            .draw()
            .unwrap()
    }
}

fn main() {
    if DO_TESTS {
        single_panel_cases();
        compound_cases();
    }
    // An alpha sweep to plot characteristics over angle of attack
    alpha_sweep();
    many_operating_points();
}
