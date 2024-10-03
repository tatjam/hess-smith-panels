extern crate plotters as pl;

use plotters::prelude::*;

pub fn plot(pts: &[(f64, f64)]) {
    let root = SVGBackend::new("out/geom.svg", (1024, 600)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .margin(15)
        .set_left_and_bottom_label_area_size(35)
        .build_cartesian_2d(-0.2..1.2, -0.1..0.25)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(std::iter::once(Polygon::new(pts, RED.mix(0.2))))
        .unwrap();

    chart
        .draw_series(
            pts.iter()
                .map(|point| Circle::new(*point, 1, BLUE.filled())),
        )
        .unwrap();

    root.present().unwrap();
}
