use crate::aero;
use crate::plots::vector_field;

const NUM_POINTS: isize = 30;
const SCALE: f64 = 0.025;

fn sample_points(start: (f64, f64), end: (f64, f64)) -> Vec<(f64, f64)> {
    assert!(start.0 < end.0);
    assert!(start.1 < end.1);

    let mut out = Vec::<(f64, f64)>::new();

    for xi in -NUM_POINTS..NUM_POINTS {
        for yi in -NUM_POINTS..NUM_POINTS {
            let x = (xi as f64) / (NUM_POINTS as f64) * (end.0 - start.0) + start.0;
            let y = (yi as f64) / (NUM_POINTS as f64) * (end.1 - start.1) + start.1;
            out.push((x, y));
        }
    }
    return out;
}

pub fn single_panel_horizontal() {
    let panel = aero::Panel {
        start: (-0.5, 0.0),
        end: (0.5, 0.0),
    };
    let range = ((-1.0, -1.0), (1.0, 1.0));
    let samples = sample_points(range.0, range.1);
    let vec_field: Vec<_> = samples
        .iter()
        .map(|point| (*point, panel.source_vel_at(*point)))
        .collect();
    vector_field::plot_vector_field("horizontal", vec_field, range, SCALE);
}

pub fn single_panel_tilted() {
    let panel = aero::Panel {
        start: (-1.0, -1.0),
        end: (1.0, 1.0),
    };
    let range = ((-2.0, -2.0), (2.0, 2.0));
    let samples = sample_points(range.0, range.1);
    let vec_field: Vec<_> = samples
        .iter()
        .map(|point| (*point, panel.source_vel_at(*point)))
        .collect();
    vector_field::plot_vector_field("tilted", vec_field, range, SCALE);
}

pub fn single_panel_backwards() {
    let panel = aero::Panel {
        start: (0.5, 0.0),
        end: (-0.5, 0.0),
    };
    let range = ((-1.0, -1.0), (1.0, 1.0));
    let samples = sample_points(range.0, range.1);
    let vec_field: Vec<_> = samples
        .iter()
        .map(|point| (*point, panel.source_vel_at(*point)))
        .collect();
    vector_field::plot_vector_field("backwards", vec_field, range, SCALE);
}
