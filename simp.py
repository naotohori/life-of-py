#!/usr/bin/env python
# coding: utf-8

# This program was originally coded by Kouichirou Nakamura as simp.py

import scipy
import scipy.linalg
import pylab
import sys

class nya(object):
    def __init__(self, coords, pdb_Bfactors):
        self.coords = coords
        self.pdb_Bfactors = pdb_Bfactors
        self.num_atoms = len(self.coords)
        
    def cmpt_graph_mats(self, cutoff=10.0):
        dist = scipy.linalg.norm
        self.adj_mat = scipy.zeros((self.num_atoms, self.num_atoms))
        self.deg_mat = [0] * self.num_atoms
        for i in range(self.num_atoms - 1):
            for j in range(i + 1, self.num_atoms):
                if dist(self.coords[i, :] - self.coords[j, :]) <= cutoff:
                    self.deg_mat[i] += 1.0
                    self.deg_mat[j] += 1.0
                    self.adj_mat[i, j] = 1.0
                    self.adj_mat[j, i] = 1.0
        self.deg_mat = scipy.diag(self.deg_mat)
        self.lap_mat = self.deg_mat - self.adj_mat
        
    def cmpt_graph_eig(self):
        self.graph_eigval, self.graph_eigvec = scipy.linalg.eigh(self.lap_mat, self.deg_mat)

    def cmpt_hessian(self):
        self.hessian = scipy.zeros((3*self.num_atoms, 3*self.num_atoms))
        for i in range(self.num_atoms - 1):
            for j in range(i + 1, self.num_atoms):
                v_ij = self.coords[j, :] - self.coords[i, :]
                d2 = sum(v_ij * v_ij)
                for a in range(3):
                    for b in range(3):
                        self.hessian[3*i + a, 3*j + b] = -v_ij[a] * v_ij[b] / d2 * self.adj_mat[i, j]
                        self.hessian[3*j + b, 3*i + a] = self.hessian[3*i + a, 3*j + b]
        for i in range(self.num_atoms):
            for a in range(3):
                for b in range(a, 3):
                    for j in range(self.num_atoms):
                        if j != i: 
                            self.hessian[3*i + a, 3*i + b] += -self.hessian[3*i + a, 3*j + b]
                    self.hessian[3*i + b, 3*i + a] = self.hessian[3*i + a, 3*i + b]

    def cmpt_en_eig(self):
        self.en_eigval, self.en_eigvec = scipy.linalg.eigh(self.hessian)

    def cmpt_inverse_hessian(self):
        self.inverse_hessian = scipy.linalg.pinv(self.hessian)

    def cmpt_Bfactors(self):
        Bfactors = [self.inverse_hessian[3*i,3*i] +
                    self.inverse_hessian[3*i+1, 3*i+1] +
                    self.inverse_hessian[3*i+2, 3*i+2] 
                    for i in range(self.num_atoms)]
        k = sum(self.pdb_Bfactors) / sum(Bfactors)
        self.Bfactors = [Bfactors[i] * k for i in range(self.num_atoms)]

    def cmpt_cross_correlation(self):
        self.cross_correlation = scipy.zeros((self.num_atoms, self.num_atoms))
        self.norm_cross_correlation = scipy.zeros((self.num_atoms, self.num_atoms))
        for i in range(self.num_atoms):
            for j in range(i, self.num_atoms):
                self.cross_correlation[i, j] = (self.inverse_hessian[3*i, 3*j] + 
                                                self.inverse_hessian[3*i+1, 3*j+1] +
                                                self.inverse_hessian[3*i+2, 3*j+2])
                self.cross_correlation[j, i] = self.cross_correlation[i, j]
        for i in range(self.num_atoms):
            for j in range(i, self.num_atoms):
                if i == j:
                    self.norm_cross_correlation[i, i] = 1.0
                else:
                    self.norm_cross_correlation[i, j] = (
                     self.cross_correlation[i, j] /
                     scipy.sqrt(self.cross_correlation[i, i] * 
                                self.cross_correlation[j, j]))
                    self.norm_cross_correlation[j, i] = self.norm_cross_correlation[i, j]
        
    
def get_lines(filename):
    lines = []
    for line in open(filename):
        if (line[0:6] == "ATOM  " and
            line[12:16] == " CA " and
            (line[16:17] == " " or line[16:17] == "A") and
            line[21:22] == "A"):
            lines.append(line)
    return lines

def get_coords(lines):
    def ext_coords(line):
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            print("Invalid format(coords)")
            print(line)
            quit()
        return (x, y, z)
    return [ext_coords(line) for line in lines]

def get_Bfactors(lines):
    def ext_Bfactors(line):
        try:
            b = float(line[60:66])
        except ValueError:
            print("Invalid format(B-factors)")
            print(line)
            quit()
        return b
    return [ext_Bfactors(line) for line in lines]

def plot_figs():
    pylab.subplot(221, aspect="equal")
    X, Y = pylab.meshgrid(list(range(ins.num_atoms)), list(range(ins.num_atoms)))
    pylab.pcolor(X, Y, ins.norm_cross_correlation)
    pylab.colorbar()
    pylab.clim(-0.15, 0.15)
    pylab.title("Cross Correlations")

    pylab.subplot(222)
    pylab.plot(pdb_Bfactors, "bo-", label="ex.")
    pylab.plot(ins.Bfactors, "ro-", label="calc.")
    pylab.legend()
    pylab.xlabel("Residue")
#    pylab.ylabel("a.u.")
    pylab.title("B factors")
    pylab.grid()

    pylab.subplot(223, aspect="equal")
    X, Y = pylab.meshgrid(list(range(ins.num_atoms)), list(range(ins.num_atoms)))
    pylab.pcolor(X, Y, ins.adj_mat)
    pylab.colorbar()
    pylab.title("Adjacency Mat.")

    pylab.subplot(224)
    pylab.plot(ins.graph_eigvec[:, 1], "go-")
    pylab.xlabel("Residue")
    pylab.grid()

    pylab.show()


if __name__ == "__main__":

    filename = sys.argv[1]

    lines = get_lines(filename)
    coords = scipy.array(get_coords(lines))
    pdb_Bfactors = get_Bfactors(lines)

    ins = nya(coords, pdb_Bfactors)
    ins.cmpt_graph_mats()
    ins.cmpt_graph_eig()
    ins.cmpt_hessian()
    ins.cmpt_en_eig()
    ins.cmpt_inverse_hessian()
    ins.cmpt_Bfactors()
    ins.cmpt_cross_correlation()

    plot_figs()

