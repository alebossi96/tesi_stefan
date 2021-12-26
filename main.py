import core.ForwardModel as fwd
import core.SignalsAndSystems as sis
import core.ResolutionToPoints as rtp
import core.Geometry as geo

geometry = geo.Geometry(thickness_top = 5, thickness_bottom = 10,
                        source_det_distance = 0,
                        mu_a_top = 0.01, mu_a_bottom = 0.01,
                        mu_s_top = 0.5, mu_s_bottom = 0.5,
                        mu_R_top = 0.011, mu_R_bottom = 0.011)
forward = fwd.ForwardModel(geometry)


instrument_data = rtp.ResolutionToPoints(resolution = 10, n_points = 1024, min_wavenumber=800, max_wavenumber = 3000)
s_top = sis.Spectrum(info = "Spectra of first layer",instrumentData =instrument_data)
s_bottom = sis.Spectrum(info = "Spectra of bottom layer",instrumentData = instrument_data)
s_top.set_gaussian(mean = 1800, amplitude = 0.012)
s_bottom.set_gaussian(mean = 1086, amplitude = 0.012)
s_top.plot(xlabel = "Raman shift $[cm^{-1}]$", ylabel ="$\mu_{s, R_1}$")
s_bottom.plot(xlabel = "Raman shift $[cm^{-1}]$", ylabel ="$\mu_{s, R_2}$")

