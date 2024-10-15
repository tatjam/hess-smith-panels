pub struct Panel {
    pub start: (f64, f64),
    pub end: (f64, f64),
}

impl Panel {
    pub fn theta(&self) -> f64 {
        //((self.end.1 - self.start.1) / (self.end.0 - self.start.0)).atan()
        //    + std::f64::consts::PI * 0.5
        (self.end.1 - self.start.1).atan2(self.end.0 - self.start.0) /*+ std::f64::consts::PI * 0.5*/
    }
    pub fn params_at(&self, p: (f64, f64)) -> (f64, f64, f64) {
        if relative_eq!(0.5 * (self.start.0 + self.end.0), p.0)
            && relative_eq!(0.5 * (self.start.1 + self.end.1), p.1)
        {
            return (0.0, std::f64::consts::PI, 0.0);
        }

        let r_ij1 = ((self.end.0 - p.0).powi(2) + (self.end.1 - p.1).powi(2)).sqrt();
        let r_ij = ((self.start.0 - p.0).powi(2) + (self.start.1 - p.1).powi(2)).sqrt();

        let betanum =
            (p.1 - self.end.1) * (p.0 - self.start.0) - (p.0 - self.end.0) * (p.1 - self.start.1);
        let betaden =
            (p.0 - self.end.0) * (p.0 - self.start.0) + (p.1 - self.end.1) * (p.1 - self.start.1);
        //let beta = (betanum / betaden).atan() + std::f64::consts::PI * 0.5;
        //let beta = betanum.atan2(betaden) + std::f64::consts::PI * 0.5;
        /*let mut beta = (betanum / betaden).atan() + std::f64::consts::PI * 0.5;
        if betanum < 0.0 {
            beta = -beta;
        }*/
        let beta = betanum.atan2(betaden);

        ((r_ij1 / r_ij).ln(), beta, 1.0)
    }

    // Induced velocity by panel on given point, assuming unitary source strength
    // (U_M,ij W_M,ij)
    pub fn source_vel_at(&self, p: (f64, f64)) -> (f64, f64) {
        let params = self.params_at(p);
        let theta_j = self.theta();

        // Panel coordinates
        let u_mij = -0.5 * std::f64::consts::FRAC_1_PI * params.0;
        let w_mij = 0.5 * std::f64::consts::FRAC_1_PI * params.1;

        // Global coordinates
        let gu_mij = u_mij * theta_j.cos() - w_mij * theta_j.sin();
        let gw_mij = u_mij * theta_j.sin() + w_mij * theta_j.cos();

        (gu_mij, gw_mij)
    }

    // Induced velocity by panel on given point, assuming unitary vortex strength
    pub fn vortex_vel_at(&self, p: (f64, f64)) -> (f64, f64) {
        let params = self.params_at(p);
        let theta_j = self.theta();

        // Panel coordinates
        let u_tij = 0.5 * std::f64::consts::FRAC_1_PI * params.1;
        let w_tij = 0.5 * std::f64::consts::FRAC_1_PI * params.0;

        // Global coordinates
        let gu_tij = u_tij * theta_j.cos() - w_tij * theta_j.sin();
        let gw_tij = u_tij * theta_j.sin() + w_tij * theta_j.cos();

        (gu_tij, gw_tij)
    }

    pub fn midpoint(&self) -> (f64, f64) {
        (
            0.5 * (self.start.0 + self.end.0),
            0.5 * (self.start.1 + self.end.1),
        )
    }

    pub fn len(&self) -> f64 {
        ((self.end.0 - self.start.0).powi(2) + (self.end.1 - self.start.1).powi(2)).sqrt()
    }

    pub fn normal(&self) -> (f64, f64) {
        let dir = ((self.end.0 - self.start.0), (self.end.1 - self.start.1));
        let dir_nrm = (dir.0 / self.len(), dir.1 / self.len());
        return (-dir_nrm.1, dir_nrm.0);
    }
}

mod test {
    #[allow(unused_imports)]
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
