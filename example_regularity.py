from classes import *
#try "bear", "dragon", "square"
tile = Tile("bear")
num_conv = 3
spline = Bspline(tile, num_conv)
print("Holder exponent of tile in L_2 is", tile.basic_smooth())
print("Basic approach for Holder exponent in L_2 gives alpha_2 =", spline.basic_smooth())
print("Spline method for Holder exponent in L_2 gives alpha_2 =", spline.smoothness())
