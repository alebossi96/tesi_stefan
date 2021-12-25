import core.ForwardModel as fwd
import core.SignalsAndSystems as sis
import core.ResolutionToPoints as rtp
instrument_data = rtp.ResolutionToPoints(resolution = 10, n_points = 1024, min_wavenumber=800, max_wavenumber = 3000)
s_top = sis.Spectrum(info = "Spectra of first layer",instrumentData =instrument_data)
s_bottom = sis.Spectrum(info = "Spectra of bottom layer",instrumentData = instrument_data)
s_top.set_gaussian(mean = 1800, amplitude = 0.012)
s_bottom.set_gaussian(mean = 1086, amplitude = 0.012)
s_top.plot(xlabel = "Raman shift $[\frac{1}{cm}]$", ylabel ="$\mu_{s, R_1}$")
s_bottom.plot(xlabel = "Raman shift $[\frac{1}{cm}]$", ylabel ="$\mu_{s, R_2}$")
forward = fwd.ForwardModel()
