from classes import *
#try "bear", "dragon", "square"
tile = Tile("bear")
num_conv = 3
spline = Bspline(tile, num_conv)
print("The Holder exponent of the tile in L_2 is", tile.basic_smooth())
print("The basic approach for Holder exponent in L_2 gives alpha_2 =", spline.basic_smooth())
print("The B-spline method for Holder exponent in L_2 gives alpha_2 =", spline.smoothness())
