extern crate plotters as pl;
use plotters::prelude::*;

pub fn plot_stream_plot<F: Fn((f64, f64)) -> (f64, f64)>(
    name: &str,
    start_points: Vec<(f64, f64)>,
    sample_vel_fn: F,
    range: ((f64, f64), (f64, f64)),
    dt: f64,
    et: f64,
    pts: Option<&[(f64, f64)]>,
) {
    let mut name_str: String = "out/".to_owned();
    name_str.push_str(name);
    name_str.push_str(".svg");
    // x / y, thus x = y * AR
    let aspect_ratio = (range.1 .0 - range.0 .0) / (range.1 .1 - range.0 .1);
    let root = SVGBackend::new(name_str.as_str(), ((512.0 * aspect_ratio) as u32, 512))
        .into_drawing_area();

    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .margin(15)
        .set_left_and_bottom_label_area_size(35)
        .build_cartesian_2d(range.0 .0..range.1 .0, range.0 .1..range.1 .1)
        .unwrap();

    // Propagate stream lines
    let streamlines: Vec<Vec<(f64, f64)>> = start_points
        .into_iter()
        .map(|(x0, y0)| {
            let mut out = Vec::<(f64, f64)>::new();
            let mut t = 0.0;
            let (mut x, mut y) = (x0, y0);
            out.push((x0, y0));
            // Simple as it can be Euler integrator
            while t < et {
                let (vx, vy) = sample_vel_fn((x, y));
                x -= vx * dt;
                y -= vy * dt;

                t += dt;
                out.push((x, y));
            }
            return out;
        })
        .collect();

    chart
        .draw_series(
            streamlines
                .into_iter()
                .map(|sline| PathElement::new(sline, BLUE.stroke_width(2))),
        )
        .unwrap();

    if pts.is_some() {
        chart
            .draw_series(std::iter::once(Polygon::new(pts.unwrap(), RED.mix(0.2))))
            .unwrap();
    }

    chart.configure_mesh().draw().unwrap();

    root.present().unwrap();
}
