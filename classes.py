import numpy as np
from math import sqrt, log
from copy import copy, deepcopy
from numpy.linalg import inv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class RefinEq():
    '''Defines refinement equation (RE)
    f(x) = c_0 f(Mx - d_0) + c_1 f(Mx - d_1) + ... + c_k f(Mx - d_k)
    Params:
    M     -- expanding integer matrix, scaling
    coef  -- list of coefficients [c_0, ..., c_k]
    digs  -- list of shift vectors [d_0, ..., d_k]
    tile  -- tile corresponding to this refinement equation (class Tile)
    listT -- transition matrices
    Omega -- minimal set such that the support of solution is a subset of tile + Omega
    N     -- len(Omega)
    reff  -- solution (class Refin_Func)'''
    def __init__(self, matrix, coefficients, coordinates):
        self.M = np.array(matrix)
        self.M1 = np.linalg.inv(matrix)
        self.m = int(abs(np.linalg.det(self.M)))
        self.coef = np.array(coefficients)
        self.digs = np.array(coordinates)
        if abs(sum(self.coef) - self.m) > 1e-10:
            print(self.coef, self.digs, self.m)
            print("WARNING: Sum of the coefficients in RE is not equal to the det of matrix")
        self.tile = None
        self.Omega = None
        self.N = None
        self.listT = None
        self.reff = None

    def add_tile(self, dig=None):
        '''automatically adds tile in simple cases'''
        if dig is None:
            if abs(self.m - 2) > 1e-8:
                print("WARNING: you probably should specify digits in add_tile")
                self.tile = Tile(self.M, np.vstack((-np.arange((self.m)), np.zeros((self.m)))).T)
            else:
                self.tile = Tile(self.M, [[0, 0], [-1, 0]])
        else:
            self.tile = Tile(self.M, dig)

    def add_Omega(self):
        '''finds minimal Omega. requires digs != tile.digs
        implementation of the algorithm by T. Meystrick'''
        if self.tile is None:
            self.add_tile()
        self.Omega = np.array([[0, 0]])
        EPS = 1e-10
        change = True
        while change:
            change = False
            now = deepcopy(self.Omega)
            for ak in self.digs:
                for d in self.tile.digs:
                    for omegai in now:
                        s = ak + omegai - d
                        r1 = self.M1.dot(s)
                        if (abs(r1[0] - round(r1[0])) < EPS and abs(r1[1] - round(r1[1])) < EPS):
                            isinomega = False
                            for u in self.Omega:
                                if abs(u[0] - r1[0]) < EPS and abs(u[1] - r1[1]) < EPS:
                                    isinomega = True
                                    break
                            if not isinomega:
                                self.Omega = np.concatenate([self.Omega, [r1]])
                                change = True
        self.N = len(self.Omega)
        if self.N == 1:
            raise ValueError("Degenerate Omega")

    def add_listT(self):
        '''finds transition matrices'''
        if self.Omega is None:
            self.add_Omega()
        self.listT = []
        for d in self.tile.digs:
            T = np.zeros((self.N, self.N))
            for i in range(self.N):
                for j in range(self.N):
                    c = self.M.dot(self.Omega[i]) - self.Omega[j] + d
                    for p in range(len(self.digs)):
                        d1 = self.digs[p]
                        if np.linalg.norm(c - d1) < 1e-10:
                            T[i][j] = self.coef[p]
            if np.linalg.norm(sum(T) - np.ones(self.N)) > 1e-8:
                raise ValueError("Sum in T != 1")
            self.listT.append(T)

    def add_solution(self, acc=8):
        '''finds solution of refinement equation with given accuracy
        in case of many digits only x,y,z-values are supported'''
        if self.Omega is None:
            self.add_Omega()
        if self.listT is None:
            self.add_listT()
        if len(self.listT) != 2:
            self.add_solution_manydigs(acc)
            print("WARNING: many-digits implementation has limitations")
            return
        blocklen = 1 << acc
        omegalen = len(self.Omega)
        leng = omegalen * blocklen
        xx = np.zeros(leng)
        yy = np.zeros(leng)
        zz = np.zeros(leng)
        ind_where_1 = np.where(abs(np.linalg.eigvals(self.listT[0])-1) < 1e-7)[0][0]
        v0 = np.linalg.eig(self.listT[0])[1].T[ind_where_1]
        v0 /= sum(v0)
        v0 /= np.sqrt(sum(v0*v0))
        dic = {}
        for j in range(omegalen):
            zz[j * blocklen] = v0[j]
            dic[(int(self.Omega[j][0]), int(self.Omega[j][1]))] = j
        for i in range(blocklen):
            prev = i // 2
            curxy = np.array([xx[prev], yy[prev]])
            xy = self.M1.dot(curxy + self.tile.digs[i % 2])
            xx[i] = xy[0]
            yy[i] = xy[1]
            vec_prev = np.zeros(omegalen)
            for j in range(omegalen):
                vec_prev[j] = zz[prev + j * blocklen]
            znow = self.listT[i % 2].dot(vec_prev)
            for j in range(omegalen):
                zz[i + j * blocklen] = znow[j]
                xx[j * blocklen + i] = xx[i] + self.Omega[j][0]
                yy[j * blocklen + i] = yy[i] + self.Omega[j][1]
        self.reff = Refin_Func(dic, acc, xx, yy, zz)

    def add_solution_manydigs(self, acc=8):
        '''utility finction. finds solution in simple form as lists with x,y,z-values'''
        if self.listT is None:
            self.add_listT()
        xx = []
        yy = []
        zz = []
        ind_where_1 = np.where(abs(np.linalg.eigvals(self.listT[0])-1) < 1e-7)[0][0]
        v0 = np.linalg.eig(self.listT[0])[1].T[ind_where_1]
        v0 /= sum(v0)
        v0 /= np.sqrt(sum(v0*v0))
        cur_z = [v0]
        cur_xy = [np.array([0, 0])]
        for i in range(acc):
            new_cur_z = []
            new_cur_xy = []
            for j in range(len(cur_z)):
                for k in range(len(self.listT)):
                    new_cur_z.append(self.listT[k].dot(cur_z[j]))
                    new_cur_xy.append(self.M1.dot(cur_xy[j] + self.tile.digs[k]))
            cur_z = deepcopy(new_cur_z)
            cur_xy = deepcopy(new_cur_xy)
        for j in range(len(self.Omega)):
            for h in range(len(cur_z)):
                xx.append(cur_xy[h][0] + self.Omega[j][0])
                yy.append(cur_xy[h][1] + self.Omega[j][1])
                zz.append(float(cur_z[h][j]))
        xx = np.array(xx)
        yy = np.array(yy)
        zz = np.array(zz)
        self.reff = Refin_Func(accur=acc, x=xx, y=yy, z=zz)

    def draw_solution(self):
        '''Plot the solution'''
        if self.reff is None:
            self.add_solution()
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(self.reff.xx, self.reff.yy, self.reff.zz, linewidth=0.2, antialiased=True)
        plt.show()

    def norm_l2(self):
        '''Normalizes solution using numerical integration'''
        if self.reff is None:
            self.add_solution()
        self.reff.zz /= np.sqrt(sum(self.reff.zz**2)/(len(self.listT) ** self.reff.acc))

    def change_signs(self):
        '''changes signs of vectors of RE'''
        return RefinEq(self.M, self.coef, -self.digs)

    def mask(self, xx, yy):
        '''finds mask of RE'''
        res = 0
        for i in range(len(self.digs)):
            k = self.digs[i][0]
            l = self.digs[i][1]
            exp = np.e ** ((-2 * np.pi * 1j) * (k * xx + l * yy))
            res += self.coef[i] * exp
        return res / self.m

    def Phi(self, xx, yy):
        '''utility periodization for construction of wavelet function'''
        sp4 = self.RE_convolution(self.change_signs())
        sp4.add_listT()
        T0 = sp4.listT[0]
        ind_where_1 = np.where(abs(np.linalg.eigvals(T0)-1) < 1e-7)[0][0]
        v0 = np.linalg.eig(T0)[1].T[ind_where_1].astype(float)
        v0 /= v0[0]
        vv = v0.astype(np.float64)
        nonz = np.where(vv != 0)
        vvnew = vv[nonz]
        Omsh = sp4.Omega[nonz]
        Omint = Omsh.astype(int)
        res = 0
        for i in range(len(Omint)):
            exp = np.e ** ((-2 * np.pi * 1j) * (Omint[i][0] * xx + Omint[i][1] * yy))
            res += vvnew[i] * exp
        return res

    def ortho_coeff_gen(self, xx, yy):
        '''utility Fourier series with ortogonalized coefficients'''
        rr = self.Phi(xx, yy)
        if len(np.where(rr < 0)[0]) > 0:
            print("WARNING: negative values when search wavelet coefficients")
        return (1 / np.sqrt(np.float64(rr)))

    def orthogonalized_mask(self, xx, yy):
        '''orthogonalizes mask'''
        Phi1 = self.Phi(xx, yy)
        xx1 = self.M.T[0][0] * xx + self.M.T[0][1] * yy
        yy1 = self.M.T[1][0] * xx + self.M.T[1][1] * yy
        Phi2 = self.Phi(xx1, yy1)
        mas = self.mask(xx, yy)
        return mas * np.sqrt(Phi1) / np.sqrt(Phi2)

    def RE_convolution(self, ref2):
        '''implements convolution of refinement equation'''
        if np.linalg.norm(self.M.ravel() - ref2.M.ravel()) > 1e-8:
            raise ValueError("Impossible to convolve two RE with different matrices")
        dic = {}
        for i in range(len(self.digs)):
            for j in range(len(ref2.digs)):
                dignew = self.digs[i] + ref2.digs[j]
                if (round(dignew[0]), round(dignew[1])) in dic:
                    dic[(round(dignew[0]), round(dignew[1]))] += self.coef[i] * ref2.coef[j]
                else:
                    dic[(round(dignew[0]), round(dignew[1]))] = self.coef[i] * ref2.coef[j]
        dig = []
        coe = []
        for el in dic:
            dig.append(np.array([el[0], el[1]]))
            coe.append(dic[el] / self.m)
        return RefinEq(self.M, coe, dig)

    def constrain_W(self):
        '''Restricts all transition matrices to the space {x_1 + ... x_k = 0}'''
        basis = np.eye(self.N)
        basis -= np.eye(self.N, self.N, 1)
        basis[-1][-1] = 0
        basis[-1][0] = 1
        list_A = []
        for T in self.listT:
            A0 = np.zeros((self.N - 1, self.N - 1))
            util = np.linalg.pinv(basis.T)
            for i in range(self.N - 1):
                vect = T.dot(basis[i])
                next_col = util.dot(vect)
                if abs(next_col[-1]) > 1e-6:
                    raise ValueError("Error: please check sum=1 in T")
                A0[i] = copy(next_col[:self.N - 1])
            A0 = A0.T
            list_A.append(A0)
        return list_A

    def calc_by_kron(self, list_A):
        '''Kronecker product's approach for joint spectral radius calculation'''
        list_kron = []
        for mA in list_A:
            list_kron.append(np.kron(mA, mA))
        list_kron = np.array(list_kron)
        final = sum(list_kron) / len(list_kron)
        other_method_ans = np.linalg.eigvals(final)
        m_other = -np.inf
        for i in range(len(other_method_ans)):
            if abs(other_method_ans[i].imag) < 1e-8:
                m_other = max(abs(other_method_ans[i].real), m_other)
        return m_other

    def basic_smooth(self):
        '''Calculates Holder exponent by
        restricting transition matrices to {x_1 + ... x_k = 0}.
        Works if there is no common invariant subspaces
        (for tiles almost always the case)'''
        if self.listT is None:
            self.add_solution()
        list_A = self.constrain_W()
        rho_sq = self.calc_by_kron(list_A)
        Meigv = np.linalg.eigvals(self.M)
        modules = (Meigv * Meigv.conjugate())**0.5
        Mrad = float(max(modules.real))
        assert(Mrad > 1)
        return round(log(sqrt(rho_sq), 1. / Mrad), 5)


class Tile(RefinEq):
    '''Defines Tile with refinement equation
    chi_G(x) = chi_G(Mx - d_0) + chi_G(Mx - d_1) + ... + chi_G(Mx - d_m), m = |det M|
    Can be defined as constructor(matrix, digits) or constructor(name),
    where name = dragon, bear or rectangle (possible tiles with 2 digits)'''
    def __init__(self, matrix, digits=-1):
        if type(digits) != int:
            super().__init__(matrix, np.ones(len(digits)), digits)
            self.check_digits()
        if type(digits) == int and digits == -1:
            if matrix == "dragon":
                super().__init__(np.array([[1, 1], [-1, 1]]), np.ones(2), np.array([[0, 0], [1, 0]]))
            if matrix == "bear":
                super().__init__(np.array([[1, -2], [1, 0]]), np.ones(2), np.array([[0, 0], [1, 0]]))
            if matrix == "square":
                super().__init__(np.array([[0, -2], [1, 0]]), np.ones(2), np.array([[0, 0], [1, 0]]))

    def check_digits(self):
        '''checks that digits in tile are from different classes Z^2 / MZ^2'''
        for dig1 in self.digs:
            for dig2 in self.digs:
                diff = dig2 - dig1
                if not(diff[0] == 0 and diff[1] == 0):
                    res = self.M1.dot(diff)
                    if abs(res[0] - round(res[0])) < 1e-10 and abs(res[1] - round(res[1])) < 1e-10:
                        raise ValueError("Bad digits in tile")
        return True


class Bspline(RefinEq):
    '''Defines B-splines as a convolution of tiles
    if typp = "symm": B-spline = chi_tile * ... * chi_tile (num_convolution)
    else: B-spline = (chi_tile * chi_{-tile})* ... * (chi_tile * chi_{-tile}) (num_convolution)'''
    def __init__(self, btile, num_convolution, typp="symm"):
        self.basetile = btile
        self.typ = typp
        self.num_conv = num_convolution
        if typp == "symm":
            spl = btile
            for i in range(self.num_conv - 1):
                spl = spl.RE_convolution(btile)
            super().__init__(spl.M, spl.coef, spl.digs)
        else:
            btile_min = Tile(btile.M, -btile.digs)
            base_elem = btile.RE_convolution(btile_min)
            spl = base_elem
            for i in range(self.num_conv - 1):
                spl = spl.RE_convolution(base_elem)
            super().__init__(spl.M, spl.coef, spl.digs)
        if len(btile.digs) > 2 and num_convolution > 1:
            self.add_tile(btile.digs)

    def generate_basis(self):
        '''generates a special basis for transition matrices'''
        nn = self.num_conv
        if self.typ != 'symm':
            nn *= 2
        basis = np.zeros((nn * (nn + 1) // 2, self.N))
        ind = 0
        for k in range(nn):
            for j in range(k + 1):
                basis[ind] = copy((self.Omega.T[0] ** (k - j)) * (self.Omega.T[1] ** j))
                ind += 1
        return basis

    def in_basis(self, vect, basold, bas, columns):
        '''represents vector from one basis in another reducing dimension'''
        nn = len(bas[0])
        short = []
        for col in columns:
            short.append(vect[col])
        bas1 = np.array(bas)
        bas1 = bas1.T
        invv = np.linalg.pinv(bas1)
        coef = invv.dot(short)
        for i in range(len(coef)):
            coef[i] = round(coef[i], 10)
            if abs(coef[i]) < 1e-10:
                coef[i] = 0
        # self-check
        our_answer = np.zeros(len(vect))
        for i in range(len(coef)):
            our_answer += coef[i] * np.array(basold[i])
        err = sum((our_answer - vect) * (our_answer - vect))
        if (err > 1e-5):
            print("WARNING: error in in_basis()", err, vect, our_answer)
        return coef

    def are_basis_with_last(self, array):
        '''checks linear independence'''
        mull = array.dot(array.T)
        if abs(np.linalg.det(mull)) < 1e-5:
            return False
        return True

    def find_orthogonal_bas(self, basis):
        '''adds vectors until basis is full'''
        num = len(basis[0])
        j = 0
        while len(basis) < num:
            vect = np.zeros(num)
            vect[j] = 1
            basis = np.vstack((basis, vect))
            if not self.are_basis_with_last(basis):
                basis = basis[:-1]
            j += 1
        return basis

    def smoothness(self):
        '''Calculates Holder exponent for B-splines'''
        basold = self.generate_basis()
        cou = int(self.num_conv * (self.num_conv + 1) * 0.5)
        bas = []
        j = 0
        columns = []
        while len(bas) < cou:
            bas.append(basold.T[j])
            if np.linalg.matrix_rank(np.array(bas)) != min(np.array(bas).shape):
                bas = bas[:-1]
            else:
                columns.append(j)
            j += 1
        bas = np.array(bas).T[:cou]
        assert(abs(np.linalg.det(bas)) > 1e-5)
        arr = []
        for vec in basold:
            arr.append(self.in_basis(self.listT[0].T.dot(vec), basold, bas, columns))
        arr = np.array(arr)
        arr = arr.T
        orth = self.find_orthogonal_bas(copy(basold))
        list_of_last_blocks = []
        for mT in self.listT:
            arr_orth = []
            for i in range(len(basold), len(orth)):
                arr_orth.append(self.in_basis(mT.T.dot(orth[i]), orth, orth, [i for i in range(len(orth))]))
            arr_orth = np.array(arr_orth)
            arr_orth = arr_orth.T
            arr_orth = arr_orth[len(basold):]
            arr_orth = arr_orth.T
            list_of_last_blocks.append(arr_orth)
        rho_2 = self.calc_by_kron(list_of_last_blocks)
        Meigv = np.linalg.eigvals(self.M)
        modules = (Meigv * Meigv.conjugate())**0.5
        Mrad = float(max(modules.real))
        assert(Mrad > 1)
        return round(log(sqrt(rho_2), 1. / Mrad), 5)


class Refin_Func():
    '''Defines refinement function on the union of several tiles.
    Values are stored in consecutive blocks with the length = blocklen.
    Shifts are defined by dictionary dic = {shift: number of block}
    The structure in block corresponds to M-nary system'''
    def __init__(self, diction=None, accur=None, x=None, y=None, z=None):
        self.dic = diction
        self.acc = accur
        self.blocklen = (1 << self.acc)
        self.xx = x
        self.yy = y
        self.zz = z

    def draw(self):
        '''Plots the values of function'''
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.view_init(elev=10, azim=130)
        ax.plot_trisurf(self.xx, self.yy, self.zz, linewidth=0.2, antialiased=True)
        plt.show()

    def shift(self, shif):
        '''shifts the function by vector shif'''
        shifted_zz = np.zeros(self.zz.shape)
        for vec in self.dic:
            vec1 = (vec[0] - shif[0], vec[1] - shif[1])
            if vec1 in self.dic:
                val = self.zz[self.dic[vec1] * self.blocklen:(self.dic[vec1] + 1) * self.blocklen]
                shifted_zz[self.dic[vec] * self.blocklen:(self.dic[vec] + 1) * self.blocklen] = val
        return Refin_Func(self.dic, self.acc, self.xx, self.yy, shifted_zz)

    def linear_combination(self, nonzero):
        '''Calculates linear combination sum c_k f(x - shift_k) defined by nonzero
        nonzero -- dictionary {shift: coeff}'''
        dic1 = self.dic.copy()
        if self.dic[0, 0] != 0:
            raise ValueError("WARNING: Omega[0] != (0, 0)")
        cur_length = len(dic1)
        newxx = self.xx
        newyy = self.yy
        newzz = np.zeros(newxx.shape)
        for (a, b) in nonzero.keys():
            a, b = int(a), int(b)
            if abs(nonzero[(a, b)]) > 1e-5:
                for om in self.dic:
                    veca = (om[0] + a, om[1] + b)
                    if veca in dic1:
                        ind = dic1[veca]
                    else:
                        dic1[veca] = cur_length
                        cur_length += 1
                        newxx = np.hstack((newxx, newxx[:self.blocklen] + veca[0]))
                        newyy = np.hstack((newyy, newyy[:self.blocklen] + veca[1]))
                        newzz = np.hstack((newzz, np.zeros(self.blocklen)))
                        ind = dic1[veca]
                    for j in range((self.blocklen)):
                        ad_val = self.zz[dic1[(om[0], om[1])] * self.blocklen + j]
                        newzz[ind * self.blocklen + j] += nonzero[(a, b)] * ad_val
        return Refin_Func(dic1, self.acc, newxx, newyy, newzz)

    def cut_zeros(self, eps=1e-2):
        '''Reduces number of points by deleting blocks with values less than eps '''
        strxx = np.array([])
        stryy = np.array([])
        strzz = np.array([])
        dic1 = {}
        ind = 0
        for (bl, i) in self.dic.items():
            if max(abs(self.zz[i * self.blocklen:(i + 1) * self.blocklen])) > eps:
                strxx = np.hstack((strxx, self.xx[i * self.blocklen:(i + 1) * self.blocklen]))
                stryy = np.hstack((stryy, self.yy[i * self.blocklen:(i + 1) * self.blocklen]))
                strzz = np.hstack((strzz, self.zz[i * self.blocklen:(i + 1) * self.blocklen]))
                dic1[bl] = ind
                ind += 1
        return Refin_Func(dic1, self.acc, strxx, stryy, strzz)

    def wavelet_func(self, nonzero1, req, formula_shift=[0, 1]):
        '''calculates wavelet function for 2-digit B-splines'''
        cur_length = len(self.dic)
        dic1 = self.dic.copy()
        wavexx = self.xx.copy()
        waveyy = self.yy.copy()
        wavezz = np.zeros(wavexx.shape)
        if dic1[0, 0] != 0:
            raise ValueError("WARNING: Omega[0] != (0, 0)")
        for (a, b) in nonzero1.keys():
            if abs(nonzero1[(a, b)]) > 1e-5:
                for om in self.dic:
                    newa = req.tile.M1.dot([om[0] + formula_shift[0] - a, om[1] + formula_shift[1] - b])
                    if abs(newa[0] - int(newa[0])) < 1e-5 and abs(newa[1] - int(newa[1])) < 1e-5:
                        veca = (int(newa[0]), int(newa[1]))
                        if veca in dic1:
                            ind = dic1[veca]
                        else:
                            dic1[veca] = cur_length
                            cur_length += 1
                            wavexx = np.hstack((wavexx, wavexx[:self.blocklen] + veca[0]))
                            waveyy = np.hstack((waveyy, waveyy[:self.blocklen] + veca[1]))
                            wavezz = np.hstack((wavezz, np.zeros(self.blocklen)))
                            ind = dic1[veca]
                        for j in range((self.blocklen + 1) // 2):
                            j1 = 2 * j
                            ad_val = self.zz[dic1[(int(om[0]), int(om[1]))] * self.blocklen + j]
                            new_ind = ind * self.blocklen + j1
                            wavezz[new_ind] += nonzero1[(a, b)] * ad_val * ((-1) ** (b - a))
                    else:
                        nwa = req.tile.M1.dot([om[0] + formula_shift[0] - a - req.tile.digs[1][0],
                                               om[1] + formula_shift[1] - b - req.tile.digs[1][1]])
                        veca = (int(nwa[0]), int(nwa[1]))
                        if veca in dic1:
                            ind = dic1[veca]
                        else:
                            dic1[veca] = cur_length
                            cur_length += 1
                            wavexx = np.hstack((wavexx, wavexx[:self.blocklen] + veca[0]))
                            waveyy = np.hstack((waveyy, waveyy[:self.blocklen] + veca[1]))
                            wavezz = np.hstack((wavezz, np.zeros(self.blocklen)))
                            ind = dic1[veca]
                        for j in range(self.blocklen):
                            j1 = 2 * j + 1
                            if j1 < self.blocklen:
                                ad_val = self.zz[dic1[(int(om[0]), int(om[1]))] * self.blocklen + j]
                                new_ind = ind * self.blocklen + j1
                                wavezz[new_ind] += nonzero1[(a, b)] * ad_val * ((-1) ** (b - a))
        return Refin_Func(dic1, self.acc, wavexx, waveyy, wavezz)

    def scalar(self, reff):
        '''calculates numerically scalar product of two refinement functions'''
        return sum(self.zz * reff.zz) / self.blocklen


def fourier_series_coeff_numpy(spline, ff, N, T=1, eps=1e-4):
    '''utility function for Fourier series'''
    t = np.linspace(0, T, 2 * N + 2, endpoint=False)
    p = np.array([[h, hh] for h in t for hh in t])
    yy = np.fft.fftn(ff(p.T[0], p.T[1]).reshape((2 * N + 2, 2 * N + 2))) / len(p)
    lenn = len(yy)
    nonzero = {}
    for i in range(yy.shape[0]):
        for j in range(yy.shape[1]):
            if abs(yy[i][j].imag) > 1e-6:
                print("WARNING: complex coefficients in Fourier transform")
            if abs(yy[i][j]) > eps:
                i2 = i
                j2 = j
                if i > (lenn - i):
                    i2 = i - lenn
                if j > (lenn - j):
                    j2 = j - lenn
                nonzero[(-i2, -j2)] = yy[i][j].real
    return nonzero
