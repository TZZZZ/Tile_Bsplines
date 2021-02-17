# Tile_Bsplines

The detailed explanation is in \textit{introduction.pdf}. 

In this repository we construct a new family of multivariate B-splines, corresponding wavelet systems, and subdivision schemes. We define B-spline as an autoconvolution of the characteristic function of a special compact set <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/058c047e1c79e7c701ffd59018a85573.svg?invert_in_darkmode" align=middle width=12.92478pt height=22.46574pt/> called ``tile''. A tile is a compact subset <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/5201385589993766eea584cd3aa6fa13.svg?invert_in_darkmode" align=middle width=12.92478pt height=22.46574pt/> of <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/435f1061aa6f25938c3c3515c083d06c.svg?invert_in_darkmode" align=middle width=18.71529pt height=27.91272pt/> possessing two properties: <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/5201385589993766eea584cd3aa6fa13.svg?invert_in_darkmode" align=middle width=12.92478pt height=22.46574pt/> is self similar by  means of several contraction affine mappings with the same linear part, i.e. <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/abbf47129927d7bb1d2deb60b4bd4d8f.svg?invert_in_darkmode" align=middle width=149.113965pt height=29.6802pt/>;  <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/5201385589993766eea584cd3aa6fa13.svg?invert_in_darkmode" align=middle width=12.92478pt height=22.46574pt/> defines a tiling of <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/435f1061aa6f25938c3c3515c083d06c.svg?invert_in_darkmode" align=middle width=18.71529pt height=27.91272pt/>, i.e., its  integer shifts cover entire space in one layer (their intersections are of Lebesgue measure zero). 
The B-spline <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/7ba90f716bfb9fd9cee3656e7e8315d0.svg?invert_in_darkmode" align=middle width=56.17623pt height=24.6576pt/> is defined as <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/727f443dc41ed7bfea6ba592443f7cc2.svg?invert_in_darkmode" align=middle width=122.2188pt height=24.6576pt/> (n convolutions), 
where <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/2c0eee5e2c85d3acc524221a8881f62b.svg?invert_in_darkmode" align=middle width=84.65985pt height=24.6576pt/> if <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/058c047e1c79e7c701ffd59018a85573.svg?invert_in_darkmode" align=middle width=42.410775pt height=22.46574pt/> and <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/fe5e14e93f9587b9fb81f2dec57a641c.svg?invert_in_darkmode" align=middle width=84.65985pt height=24.6576pt/> otherwise. 

The most popular and the most practically useful are 2-tiles, which are defined by two affine contractions. On the two-dimensional plane \mathbb{R}^2,  there are exactly three tiles, up to affine similarity. They are  called  "square" (this is a real square), "dragon", and "bear" (also known as "twindragon" or  "tame twindragon"). We call corresponding B-splines Bear-k (convolution of k Bears) and similarly with the square and with the dragon. The 2-tiles can be constructed in arbitrary dimension. The Square-k is the classical B-spline, which is a direct product of <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/2103f85b8b1477f430fc407cad462224.svg?invert_in_darkmode" align=middle width=8.556075pt height=22.83138pt/> univariate splines. This is a piecewise-polynomial function, but Dragon-k and  Bear-k are neither piecewise-polynomial nor piecewise analytic. 

In the case of 2-tiles the number of coefficients of the subdivision scheme defined by each of these three B-splines of order k is <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/6b44835ef9c9df90c1ab13fe002f5bf9.svg?invert_in_darkmode" align=middle width=37.385865pt height=22.83138pt/>. In particular, the subdivision scheme defined by the classical B-spline (Square-k) also has <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/ffe84ac085f6c35a98ccc33c363f07a3.svg?invert_in_darkmode" align=middle width=50.17122pt height=24.6576pt/> coefficients. 
On the other hand, if we consider that scheme as a direct product of <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/33359de825e43daa97171e27f6558ae9.svg?invert_in_darkmode" align=middle width=37.385865pt height=22.83138pt/> univariate B-spline schemes, then the number of coefficients is <img src="https://cdn.jsdelivr.net/gh/TZZZZ/Tile_Bsplines@main/svgs/a26be85c27a95007a180dd5609fd2030.svg?invert_in_darkmode" align=middle width=57.014265pt height=27.91272pt/> (with the same generated surface!). This  is one of advantages of using B-splines defined as convolutions of two-digit tiles. 
Basically, B-splines defined this way are much simpler and more efficient in both subdivision schemes and wavelets. 

In addition, it turns out that Bear-k have higher regularities than the corresponding 
classical B-splines of the same order. For example, Bear-3 is a C^2 functions while Square-3 is not; Bear-4 is in C^3, white  Square-4 is not.  This implies that the corresponding subdivision schemes have higher rate of convergence. This also speeds up the algorithms of wavelet decompositions, of computing the wavelet coefficients, etc.  

The structure of the repository:
- classes.py -- the main file with classes: Tile, B-spline, refinement equation and their methods.  
- example_wavelet_func.py -- how to obtain B-spline, orthogonalize it, and construct wavelet function. 
- example_regularity.py -- how to compute Holder regularity of B-splines in C and in L_2. 
The code was tested only for 2-tiles. 

More about tiles, multivariate subdivision algorithms and Haar systems: 

- J. Lagarias, Y. Wang, Integral self-affine tiles in R^n. II. Lattice tilings, J. Fourier
Anal. Appl. 3 (1997), no. 1, 83 -- 102. 

- K. Grochenig, W.R. Madych, Multiresolution analysis, Haar bases, and self-similar
tilings of R^n, IEEE Trans. Inform. Theory 38 (1992), no. 2, 556 -- 568. 

- K. Grochenig, A. Haas, Self-similar lattice tilings, J. Fourier Anal. Appl., 1 (1994), no.
2, 131 -- 170. 



More about Holder regularity: 

- M. Charina, V.Yu. Protasov, Regularity of anisotropic refinable functions, Appl. Comput.
Harmon. Anal., 47 (2019), no. 3, 795 -- 821. 

- V.Yu. Protasov, The generalized spectral radius. A geometric approach, Izvestiya
Math., 61 (1997), 995 -- 1030.

