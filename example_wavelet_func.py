from classes import *
tile = Tile("bear")
initial_spline = Bspline(tile, 2)
initial_spline.add_solution(acc=7)
initial_spline.norm_l2()
# plot the initial B-spline
initial_spline.draw_solution()
orthogonal_spl_coeff = fourier_series_coeff_numpy(initial_spline,
                                                  initial_spline.ortho_coeff_gen, 50)
wavelet_coeff = fourier_series_coeff_numpy(initial_spline, initial_spline.orthogonalized_mask, 50)
for elem in wavelet_coeff:
    wavelet_coeff[elem] *= 2
ortho_spline_full = initial_spline.reff.linear_combination(orthogonal_spl_coeff)
assert(abs(ortho_spline_full.scalar(ortho_spline_full) - 1) < 5 * 1e-2)
assert(abs(ortho_spline_full.scalar(ortho_spline_full.shift([1, 2]))) < 5 * 1e-2)
ortho_spline = ortho_spline_full.cut_zeros()
# plot the ortogonalized B-spline
ortho_spline.draw()
# for Square spline it is better to use formula_shift = [1, 0]
wavelet_long = ortho_spline.wavelet_func(wavelet_coeff, initial_spline, formula_shift=[0, 1])
wavelet = wavelet_long.cut_zeros()
assert(abs(wavelet.scalar(wavelet) - 1) < 5 * 1e-2)
assert(abs(wavelet.scalar(wavelet.shift([1, 0]))) < 5 * 1e-2)
# plot the wavelet function
wavelet.draw()
