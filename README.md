# Tile_Bsplines

The detailed explanation is in introduction.pdf. 

In this repository we construct a new family of multivariate B-splines, corresponding wavelet systems, and subdivision schemes. We define B-spline as an autoconvolution of the characteristic function of a special compact set ![](https://latex.codecogs.com/gif.latex?G) called ``tile''. A tile is a compact subset ![](https://latex.codecogs.com/gif.latex?G) of $R^d$ possessing two properties: ![](https://latex.codecogs.com/gif.latex?G) is self similar by  means of several contraction affine mappings with the same linear part, i.e. ![](https://latex.codecogs.com/gif.latex?G&space;=&space;\bigcup&space;\limits_{k}&space;M^{-1}(G&space;&plus;&space;d_k);  $G$ defines a tiling of ![](https://latex.codecogs.com/gif.latex?\mathbb{R}^d), i.e., its  integer shifts cover entire space in one layer (their intersections are of Lebesgue measure zero). 
The B-spline $B(G, n)$ is defined as ![](https://latex.codecogs.com/gif.latex?\chi(G)&space;*&space;\cdots&space;*&space;\chi(G)), 
where $\chi(G)[x] = 1$ if $x \in G$ and $\chi(G)[x] = 0$ otherwise. 

The most popular and the most practically useful are 2-tiles, which are defined by two affine contractions. On the two-dimensional plane \mathbb{R}^2,  there are exactly three tiles, up to affine similarity. They are  called  "square" (this is a real square), "dragon", and "bear" (also known as "twindragon" or  "tame twindragon"). We call corresponding B-splines Bear-k (convolution of k Bears) and similarly with the square and with the dragon. The 2-tiles can be constructed in arbitrary dimension. The Square-k is the classical B-spline, which is a direct product of $d$ univariate splines. This is a piecewise-polynomial function, but Dragon-k and  Bear-k are neither piecewise-polynomial nor piecewise analytic. 

In the case of 2-tiles the number of coefficients of the subdivision scheme defined by each of these three B-splines of order k is $k + 1$. In particular, the subdivision scheme defined by the classical B-spline (Square-k) also has $(k + 1)$ coefficients. 
On the other hand, if we consider that scheme as a direct product of $k+1$ univariate B-spline schemes, then the number of coefficients is $(k + 1)^d$ (with the same generated surface!). This  is one of advantages of using B-splines defined as convolutions of two-digit tiles. 
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

