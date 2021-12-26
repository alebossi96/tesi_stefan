class Geometry:
    def __init__(self, thickness_top, thickness_bottom,
                        source_det_distance,
                        mu_a_top, mu_a_bottom,
                        mu_s_top, mu_s_bottom,
                        mu_R_top, mu_R_bottom):
        self.thickness_top = thickness_top
        self.thickness_bottom = thickness_bottom
        self.source_det_distance = source_det_distance
        self.mu_a_top = mu_a_top
        self.mu_a_bottom = mu_a_bottom
        self.mu_s_top = mu_s_top
        self.mu_s_bottom = mu_s_bottom
        self.mu_R_top = mu_R_top
        self.mu_R_bottom = mu_R_bottom
        
