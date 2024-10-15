use crate::aero;
use crate::plots::vector_field;
use crate::sample_points;

const SCALE: f64 = 0.05;

pub fn single_panel_horizontal() {
    let panel = aero::Panel {
        start: (-0.5, 0.0),
        end: (0.5, 0.0),
    };
    let range = ((-1.0, -1.0), (1.0, 1.0));
    let samples = sample_points::sample_points(range.0, range.1);
    let vec_field: Vec<_> = samples
        .iter()
        .map(|point| (*point, panel.source_vel_at(*point)))
        .collect();
    vector_field::plot_vector_field("horizontal", vec_field, range, SCALE);
}

pub fn single_panel_tilted() {
    let panel = aero::Panel {
        start: (-0.5, -0.5),
        end: (0.5, 0.5),
    };
    let range = ((-1.0, -1.0), (1.0, 1.0));
    let samples = sample_points::sample_points(range.0, range.1);
    let vec_field: Vec<_> = samples
        .iter()
        .map(|point| (*point, panel.source_vel_at(*point)))
        .collect();
    vector_field::plot_vector_field("tilted", vec_field, range, SCALE);
}

pub fn single_panel_tilted_backwards() {
    let panel = aero::Panel {
        start: (0.5, 0.5),
        end: (-0.5, -0.5),
    };
    let range = ((-1.0, -1.0), (1.0, 1.0));
    let samples = sample_points::sample_points(range.0, range.1);
    let vec_field: Vec<_> = samples
        .iter()
        .map(|point| (*point, panel.source_vel_at(*point)))
        .collect();
    vector_field::plot_vector_field("tilted_backwards", vec_field, range, SCALE);
}

pub fn single_panel_backwards() {
    let panel = aero::Panel {
        start: (0.5, 0.0),
        end: (-0.5, 0.0),
    };
    let range = ((-1.0, -1.0), (1.0, 1.0));
    let samples = sample_points::sample_points(range.0, range.1);
    let vec_field: Vec<_> = samples
        .iter()
        .map(|point| (*point, panel.source_vel_at(*point)))
        .collect();
    vector_field::plot_vector_field("backwards", vec_field, range, SCALE);
}

pub fn single_panel_vertical() {
    let panel = aero::Panel {
        start: (0.0, -0.5),
        end: (0.0, 0.5),
    };
    let range = ((-1.0, -1.0), (1.0, 1.0));
    let samples = sample_points::sample_points(range.0, range.1);
    let vec_field: Vec<_> = samples
        .iter()
        .map(|point| (*point, panel.source_vel_at(*point)))
        .collect();
    vector_field::plot_vector_field("vertical", vec_field, range, SCALE);
}

pub fn single_panel_vertical_backwards() {
    let panel = aero::Panel {
        start: (0.0, 0.5),
        end: (0.0, -0.5),
    };
    let range = ((-1.0, -1.0), (1.0, 1.0));
    let samples = sample_points::sample_points(range.0, range.1);
    let vec_field: Vec<_> = samples
        .iter()
        .map(|point| (*point, panel.source_vel_at(*point)))
        .collect();
    vector_field::plot_vector_field("vertical_backwards", vec_field, range, SCALE);
}
