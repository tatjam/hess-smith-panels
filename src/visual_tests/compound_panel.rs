use crate::aero;
use crate::plots::vector_field;
use crate::sample_points;

const SCALE: f64 = 0.025;

pub fn compound_rhombus() {
    let mut panels = Vec::<aero::Panel>::new();

    panels.push(aero::Panel {
        start: (-0.4, 0.0),
        end: (0.0, 0.4),
    });
    panels.push(aero::Panel {
        start: (0.0, 0.4),
        end: (0.4, 0.0),
    });
    panels.push(aero::Panel {
        start: (0.4, 0.0),
        end: (0.0, -0.4),
    });
    panels.push(aero::Panel {
        start: (0.0, -0.4),
        end: (-0.4, 0.0),
    });

    let range = ((-1.0, -1.0), (1.0, 1.0));
    let samples = sample_points::sample_points(range.0, range.1);
    let vec_field: Vec<_> = samples
        .iter()
        .map(|point| {
            let mut sum = (0.0, 0.0);
            for panel in &panels {
                let contrib = panel.source_vel_at(*point);
                sum.0 += contrib.0;
                sum.1 += contrib.1;
            }
            (*point, sum)
        })
        .collect();
    vector_field::plot_vector_field("rhombus", vec_field, range, SCALE);
}

pub fn compound_square() {
    let mut panels = Vec::<aero::Panel>::new();

    panels.push(aero::Panel {
        start: (-0.4, -0.4),
        end: (-0.4, 0.4),
    });
    panels.push(aero::Panel {
        start: (-0.4, 0.4),
        end: (0.4, 0.4),
    });
    panels.push(aero::Panel {
        start: (0.4, 0.4),
        end: (0.4, -0.4),
    });
    panels.push(aero::Panel {
        start: (0.4, -0.4),
        end: (-0.4, -0.4),
    });

    let range = ((-1.0, -1.0), (1.0, 1.0));
    let samples = sample_points::sample_points(range.0, range.1);
    let vec_field: Vec<_> = samples
        .iter()
        .map(|point| {
            let mut sum = (0.0, 0.0);
            for panel in &panels {
                let contrib = panel.source_vel_at(*point);
                sum.0 += contrib.0;
                sum.1 += contrib.1;
            }
            (*point, sum)
        })
        .collect();
    vector_field::plot_vector_field("square", vec_field, range, SCALE);
}
