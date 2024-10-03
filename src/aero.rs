use fa::modules::core::ComplexField;

pub struct Panel {
    start: (f64, f64),
    end: (f64, f64),
}

impl Panel {
    // Induced velocity by panel on given point, assuming unitary source strength
    pub fn source_vel_at(&self, p: (f64, f64)) -> (f64, f64) {
        let rstart = ((self.start.0 - p.0).powi(2) + (self.start.1 - p.1).powi(2)).sqrt();
        let rend = ((self.end.0 - p.0).powi(2) + (self.end.1 - p.1).powi(2)).sqrt();
        let len =
            ((self.end.0 - self.start.0).powi(2) + (self.end.1 - self.start.1).powi(2)).sqrt();

        let lx = self.end.0 - self.start.0;
        let ly = self.end.1 - self.start.1;

        let sin_theta = lx / ((1.0 + (lx * lx) / (ly * ly)).sqrt() * ly);
        let cos_theta = 1.0 / ((1.0 + (lx * lx) / (ly * ly)).sqrt());
        let ty = -p.0 * sin_theta + p.1 * cos_theta;
        if approx::relative_eq!(ty, 0.0) {
            return (0.0, 0.0);
        }

        let vel_parallel = (rstart / rend).sqrt().ln();
        // TODO: Check if atan2 should be used
        let vel_perp = -(len / ty).atan();

        // Transform back
        (
            vel_parallel * cos_theta + vel_perp * sin_theta,
            vel_perp * cos_theta - vel_parallel * sin_theta,
        )
    }

    // Induced velocity by panel on given point, assuming unitary vortex strength
    pub fn vortex_vel_at(&self, p: (f64, f64)) -> (f64, f64) {
        let as_source = self.source_vel_at(p);
        (-as_source.1, as_source.0)
    }

    pub fn midpoint(&self) -> (f64, f64) {
        (
            0.5 * (self.start.0 + self.end.0),
            0.5 * (self.start.1 + self.end.1),
        )
    }
}

mod test {
    use super::*;

    #[test]
    fn test_panel_self_influence_null() {
        let panel = Panel {
            start: (0.0, 0.0),
            end: (1.0, 0.0),
        };

        let vcenter = panel.source_vel_at((0.5, 0.0));
        approx::assert_relative_eq!(vcenter.0, 0.0);
        approx::assert_relative_eq!(vcenter.1, 0.0);
    }
}
