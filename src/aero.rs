use fa::modules::core::ComplexField;

pub struct Panel {
    pub start: (f64, f64),
    pub end: (f64, f64),
}

impl Panel {
    // Induced velocity by panel on given point, assuming unitary source strength
    pub fn source_vel_at(&self, p: (f64, f64)) -> (f64, f64) {
        // Avoid singularity
        if approx::relative_eq!(p.0, self.start.0) && approx::relative_eq!(p.1, self.start.1) {
            return (0.0, 0.0);
        }
        if approx::relative_eq!(p.0, self.end.0) && approx::relative_eq!(p.1, self.end.1) {
            return (0.0, 0.0);
        }

        let rstart = ((self.start.0 - p.0).powi(2) + (self.start.1 - p.1).powi(2)).sqrt();
        let rend = ((self.end.0 - p.0).powi(2) + (self.end.1 - p.1).powi(2)).sqrt();
        let len =
            ((self.end.0 - self.start.0).powi(2) + (self.end.1 - self.start.1).powi(2)).sqrt();

        let lx = self.end.0 - self.start.0;
        let ly = self.end.1 - self.start.1;

        // sin(arctan of (ly / lx)) approaches +-1 with lx approaching 0
        let sin_theta = if lx == 0.0 {
            if ly > 0.0 {
                1.0
            } else {
                -1.0
            }
        } else {
            ly / ((1.0 + (ly * ly) / (lx * lx)).sqrt() * lx)
        };
        let cos_theta = 1.0 / ((1.0 + (ly * ly) / (lx * lx)).sqrt());
        let ty = -(p.0 - self.start.0) * sin_theta + (p.1 - self.start.1) * cos_theta;

        let vel_parallel = (rstart / rend).sqrt().ln();
        let vel_perp = if approx::relative_eq!(ty, 0.0) {
            0.0
        } else {
            // atan2 is incorrect in this case
            -(len / ty).atan()
        };

        // Transform back
        (
            vel_parallel * cos_theta + vel_perp * sin_theta,
            -vel_perp * cos_theta + vel_parallel * sin_theta,
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
    fn panel_self_influence_null() {
        let panel = Panel {
            start: (0.0, 0.0),
            end: (1.0, 0.0),
        };

        let vcenter = panel.source_vel_at((0.5, 0.0));
        approx::assert_relative_eq!(vcenter.0, 0.0);
        approx::assert_relative_eq!(vcenter.1, 0.0);
    }

    #[test]
    fn panel_self_influence_null_tilted() {
        let panel = Panel {
            start: (1.5, 1.5),
            end: (2.5, 2.5),
        };

        let vcenter = panel.source_vel_at((2.0, 2.0));
        approx::assert_relative_eq!(vcenter.0, 0.0);
        approx::assert_relative_eq!(vcenter.1, 0.0);
    }

    #[test]
    fn panel_self_influence_null_vertical() {
        let panel = Panel {
            start: (0.0, 0.0),
            end: (0.0, 1.0),
        };

        let vcenter = panel.source_vel_at((0.0, 0.5));
        approx::assert_relative_eq!(vcenter.0, 0.0);
        approx::assert_relative_eq!(vcenter.1, 0.0);
    }

    #[test]
    fn panel_logical() {
        let panel = Panel {
            start: (0.0, 0.0),
            end: (1.0, 0.0),
        };

        let vtop = panel.source_vel_at((0.5, 1.0));
        approx::assert_relative_eq!(vtop.0, 0.0);
        assert!(vtop.1 > 0.0);

        let vbottom = panel.source_vel_at((0.5, -1.0));
        approx::assert_relative_eq!(vtop.0, 0.0);
        assert!(vbottom.1 < 0.0);

        let vtopright = panel.source_vel_at((1.0, 1.0));
        assert!(vtopright.0 > 0.0);
        assert!(vtopright.1 > 0.0);

        let vtopleft = panel.source_vel_at((0.0, 1.0));
        assert!(vtopleft.0 < 0.0);
        assert!(vtopleft.1 > 0.0);

        let vbottomleft = panel.source_vel_at((0.0, -1.0));
        assert!(vbottomleft.0 < 0.0);
        assert!(vbottomleft.1 < 0.0);

        let vbottomright = panel.source_vel_at((1.0, -1.0));
        assert!(vbottomright.0 > 0.0);
        assert!(vbottomright.1 < 0.0);

        let vleft = panel.source_vel_at((-1.0, 0.0));
        approx::assert_relative_eq!(vleft.1, 0.0);
        assert!(vleft.0 < 0.0);

        let vright = panel.source_vel_at((2.0, 0.0));
        approx::assert_relative_eq!(vright.1, 0.0);
        assert!(vright.0 > 0.0);
    }

    #[test]
    fn panel_logical_vertical() {
        let panel = Panel {
            start: (0.0, 0.0),
            end: (0.0, 1.0),
        };

        let vleft = panel.source_vel_at((-1.0, 0.5));
        approx::assert_relative_eq!(vleft.1, 0.0);
        assert!(vleft.0 < 0.0);

        let vright = panel.source_vel_at((1.0, 0.5));
        approx::assert_relative_eq!(vright.1, 0.0);
        assert!(vright.0 > 0.0);

        let vtopright = panel.source_vel_at((1.0, 1.0));
        assert!(vtopright.0 > 0.0);
        assert!(vtopright.1 > 0.0);

        let vtopleft = panel.source_vel_at((-1.0, 1.0));
        assert!(vtopleft.0 < 0.0);
        assert!(vtopleft.1 > 0.0);

        let vbottomleft = panel.source_vel_at((-1.0, -1.0));
        assert!(vbottomleft.0 < 0.0);
        assert!(vbottomleft.1 < 0.0);

        let vbottomright = panel.source_vel_at((1.0, -1.0));
        assert!(vbottomright.0 > 0.0);
        assert!(vbottomright.1 < 0.0);

        let vtop = panel.source_vel_at((0.0, 2.0));
        approx::assert_relative_eq!(vtop.0, 0.0);
        assert!(vtop.1 > 0.0);

        let vbottom = panel.source_vel_at((0.0, -1.0));
        approx::assert_relative_eq!(vbottom.0, 0.0);
        assert!(vbottom.1 < 0.0);
    }

    #[test]
    fn panel_logical_reversed() {
        let panel = Panel {
            start: (1.0, 0.0),
            end: (0.0, 0.0),
        };

        let vtop = panel.source_vel_at((0.5, 1.0));
        approx::assert_relative_eq!(vtop.0, 0.0);
        assert!(vtop.1 > 0.0);

        let vbottom = panel.source_vel_at((0.5, -1.0));
        approx::assert_relative_eq!(vtop.0, 0.0);
        assert!(vbottom.1 < 0.0);

        let vtopright = panel.source_vel_at((1.0, 1.0));
        assert!(vtopright.0 > 0.0);
        assert!(vtopright.1 > 0.0);

        let vtopleft = panel.source_vel_at((0.0, 1.0));
        assert!(vtopleft.0 < 0.0);
        assert!(vtopleft.1 > 0.0);

        let vbottomleft = panel.source_vel_at((0.0, -1.0));
        assert!(vbottomleft.0 < 0.0);
        assert!(vbottomleft.1 < 0.0);

        let vbottomright = panel.source_vel_at((1.0, -1.0));
        assert!(vbottomright.0 > 0.0);
        assert!(vbottomright.1 < 0.0);

        let vleft = panel.source_vel_at((-1.0, 0.0));
        approx::assert_relative_eq!(vleft.1, 0.0);
        assert!(vleft.0 < 0.0);

        let vright = panel.source_vel_at((2.0, 0.0));
        approx::assert_relative_eq!(vright.1, 0.0);
        assert!(vright.0 > 0.0);
    }
}
