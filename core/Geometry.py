import numpy as np

class Geometry:
    def __init__(self, thickness_top, thickness_bottom,
                        mu_a_top, mu_a_bottom,
                        mu_s_top, mu_s_bottom,
                        mu_R_top, mu_R_bottom,
                        r0=0, r1 = 3, n_step_r=5,
                        source_det_distance = 0):
        self.thickness_top = thickness_top
        self.thickness_bottom = thickness_bottom
        self.source_det_distance = source_det_distance
        self.mu_a_top = mu_a_top
        self.mu_a_bottom = mu_a_bottom
        self.mu_s_top = mu_s_top
        self.mu_s_bottom = mu_s_bottom
        self.mu_R_top = mu_R_top
        self.mu_R_bottom = mu_R_bottom
        self.r0 = r0
        self.r1=r1 #mm
        self.n_step_r = n_step_r
        (self.r, self.dr) = np.linspace(self.r0,self.r1, self.n_step_r, retstep = True)
    def __iter__(self):
        return iter(self.r.tolist())
        
