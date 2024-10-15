extern crate plotters as pl;
use plotters::prelude::*;

pub fn plot_vector_field(
    name: &str,
    vecs: Vec<((f64, f64), (f64, f64))>,
    range: ((f64, f64), (f64, f64)),
    scale: f64,
) {
    let mut name_str: String = "out/".to_owned();
    name_str.push_str(name);
    name_str.push_str(".png");
    // x / y, thus x = y * AR
    let aspect_ratio = (range.1 .0 - range.0 .0) / (range.1 .1 - range.0 .1);
    let root = BitMapBackend::new(name_str.as_str(), ((1024.0 * aspect_ratio) as u32, 1024))
        .into_drawing_area();

    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .margin(15)
        .set_left_and_bottom_label_area_size(35)
        .build_cartesian_2d(range.0 .0..range.1 .0, range.0 .1..range.1 .1)
        .unwrap();

    chart
        .draw_series(vecs.iter().map(|((x0, y0), (dx, dy))| {
            PathElement::new(
                vec![(*x0, *y0), (*x0 + (*dx * scale), *y0 + (*dy * scale))],
                BLUE.stroke_width(2),
            )
        }))
        .unwrap();

    chart
        .draw_series(
            vecs.iter()
                .map(|(pos, _)| Circle::new(*pos, 2, BLUE.filled())),
        )
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    root.present().unwrap();
}
