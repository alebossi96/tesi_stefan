class ResolutionToPoints:
    def __init__(self, resolution, n_points, min_wavenumber, max_wavenumber):
        if (n_points <= 0 or not(n_points==(int)(n_points))):
            raise(ValueError,"Positive integer expected for number of points.")
        self.sigma = resolution*n_points/(max_wavenumber - min_wavenumber)
        self.n_points = n_points
        self.min_wavenumber = min_wavenumber
        self.max_wavenumber = max_wavenumber
        
