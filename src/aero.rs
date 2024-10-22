use fa::prelude::*;
use fa::FaerMat;
use fa::Mat;
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
            return (0.0, std::f64::consts::PI, 1.0);
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
        let gu_tij = (u_tij * theta_j.cos() - w_tij * theta_j.sin()) * params.2;
        let gw_tij = (u_tij * theta_j.sin() + w_tij * theta_j.cos()) * params.2;

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

pub fn hess_smith(panels: &[Panel], u_infty: f64, alpha: f64) -> fa::Mat<f64> {
    let mut mat = fa::Mat::<f64>::zeros(panels.len() + 1, panels.len() + 1);
    let mut rhs = fa::Mat::<f64>::zeros(panels.len() + 1, 1);
    for (effect_idx, effect) in panels.iter().enumerate() {
        let theta_i = effect.theta();
        let midpoint = effect.midpoint();

        for (cause_idx, cause) in panels.iter().enumerate() {
            let theta_j = cause.theta();
            let params = cause.params_at(midpoint);
            // Source is scaled by sigma_cause (one for each panel)
            mat[(effect_idx, cause_idx)] = 0.5
                * std::f64::consts::FRAC_1_PI
                * (params.0 * (theta_i - theta_j).sin() + params.1 * (theta_i - theta_j).cos());
            // Vortex is scaled by "global" vortex intensity
            // (cause panels.len() represents the vortex)
            mat[(effect_idx, panels.len())] += 0.5
                * std::f64::consts::FRAC_1_PI
                * (-params.1 * (theta_i - theta_j).sin() + params.0 * (theta_i - theta_j).cos())
                * params.2;
        }

        // Right hand side
        rhs[(effect_idx, 0)] = -u_infty * (alpha - theta_i).sin();
    }

    // Kutta condition
    // TODO: Check this or the other way around
    let top = panels.last().unwrap();
    let bottom = panels.first().unwrap();
    let theta_n = top.theta();
    let theta_1 = bottom.theta();
    for (cause_idx, cause) in panels.iter().enumerate() {
        let params_n = cause.params_at(top.midpoint());
        let params_1 = cause.params_at(bottom.midpoint());
        let theta_j = cause.theta();

        mat[(panels.len(), cause_idx)] = 0.5
            * std::f64::consts::FRAC_1_PI
            * (params_1.1 * (theta_1 - theta_j).sin() - params_1.0 * (theta_1 - theta_j).cos()
                + params_n.1 * (theta_n - theta_j).sin()
                - params_n.0 * (theta_n - theta_j).cos());

        mat[(panels.len(), panels.len())] += 0.5
            * std::f64::consts::FRAC_1_PI
            * (params_1.1 * (theta_1 - theta_j).cos() + params_1.0 * (theta_1 - theta_j).sin())
            * params_1.2;
        mat[(panels.len(), panels.len())] += 0.5
            * std::f64::consts::FRAC_1_PI
            * (params_n.1 * (theta_n - theta_j).cos() + params_n.0 * (theta_n - theta_j).sin())
            * params_n.2;
    }

    rhs[(panels.len(), 0)] = -u_infty * ((alpha - theta_1).cos() + (alpha - theta_n).cos());

    let plu = mat.full_piv_lu();
    let x = plu.solve(&rhs);

    return x;
}

pub fn get_vel_at(
    at: (f64, f64),
    sol: &faer::Mat<f64>,
    panels: &[Panel],
    u_infty: f64,
    alpha: f64,
) -> (f64, f64) {
    let mut vel = (alpha.cos() * u_infty, alpha.sin() * u_infty);
    //let mut vel = (0.0, 0.0);
    let vortex_strength = sol[(panels.len(), 0)];
    for (idx, panel) in panels.iter().enumerate() {
        let source_strength = sol[(idx, 0)];

        let svel = panel.source_vel_at(at);
        let tvel = panel.vortex_vel_at(at);

        vel.0 += source_strength * svel.0 + vortex_strength * tvel.0;
        vel.1 += source_strength * svel.1 + vortex_strength * tvel.1;
    }

    vel
}

// Returns total cL, cDi and cM around given point
pub fn find_coeffs(
    panels: &[Panel],
    x: &fa::Mat<f64>,
    u_infty: f64,
    alpha: f64,
    moment_point: (f64, f64),
) -> (f64, f64, f64) {
    let mut cl = 0.0;
    let mut cdi = 0.0;
    let mut cm = 0.0;

    for p in panels {
        let params = p.params_at((0.0, 0.1));
        let cvel = get_vel_at(p.midpoint(), &x, &panels, u_infty, alpha);
        let nrm0 = cvel.0 * p.normal().0;
        let nrm1 = cvel.1 * p.normal().1;
        assert_relative_eq!(nrm0, -nrm1, max_relative = 0.001);

        let tvelmag2 = cvel.0 * cvel.0 + cvel.1 * cvel.1;

        let cpi = 1.0 - tvelmag2 / u_infty.powi(2);

        // We don't divide by chord, but it's 1 in the study profile!
        cl -= cpi * p.len() * (p.theta() - alpha).cos();
        cdi -= cpi * p.len() * (p.theta() - alpha).sin();
        let panel_to_moment = (
            moment_point.0 - p.midpoint().0,
            moment_point.1 - p.midpoint().1,
        );
        let panel_to_moment_len = (panel_to_moment.0.powi(2) + panel_to_moment.1.powi(2)).sqrt();
        let delta_i = panel_to_moment.1.atan2(panel_to_moment.0);
        cm -= cpi * p.len() * panel_to_moment_len * (delta_i - p.theta()).cos();
    }

    return (cl, cdi, cm);
}

pub fn vel_field_panel(
    sol: &faer::Mat<f64>,
    panels: &[Panel],
    v_infty: f64,
    alpha: f64,
) -> Vec<((f64, f64), (f64, f64))> {
    let points = panels.iter().map(|panel| panel.midpoint());

    points
        .into_iter()
        .filter_map(|point| {
            let mut vel = get_vel_at(point, sol, panels, v_infty, alpha);
            if vel.0 * vel.0 + vel.1 * vel.1 > 5.0 {
                None
            } else {
                Some((point, vel))
            }
        })
        .collect()
}

pub fn vel_field_solved(
    sol: &faer::Mat<f64>,
    panels: &[Panel],
    u_infty: f64,
    alpha: f64,
) -> Vec<((f64, f64), (f64, f64))> {
    let points = super::sample_points((-0.2, -0.1), (1.2, 0.25));

    points
        .into_iter()
        .filter_map(|point| {
            let mut vel = get_vel_at(point, sol, panels, u_infty, alpha);
            if vel.0 * vel.0 + vel.1 * vel.1 > 5.0 {
                None
            } else {
                Some((point, vel))
            }
        })
        .collect()
}

mod test {
    #[allow(unused_imports)]
    use super::*;

    #[test]
    fn panel_self_influence() {
        let panel = Panel {
            start: (0.0, 0.0),
            end: (1.0, 0.0),
        };

        let vcenter = panel.source_vel_at((0.5, 0.0));
        approx::assert_relative_eq!(vcenter.0, 0.0);
        approx::assert_relative_eq!(vcenter.1, 0.5);
    }

    #[test]
    fn panel_self_influence_tilted() {
        let panel = Panel {
            start: (1.5, 1.5),
            end: (2.5, 2.5),
        };

        let vcenter = panel.source_vel_at((2.0, 2.0));
        approx::assert_relative_eq!(vcenter.0, -vcenter.1);
    }

    #[test]
    fn panel_self_influence_vertical() {
        let panel = Panel {
            start: (0.0, 0.0),
            end: (0.0, 1.0),
        };

        let vcenter = panel.source_vel_at((0.0, 0.5));
        approx::assert_relative_eq!(vcenter.0, -0.5);
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
