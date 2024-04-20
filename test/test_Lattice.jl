using OMEinsum
using ClassicalSpinMonteCarlo
using LinearAlgebra

a1 = (1, 0)
a2 = (0, 1)
uc = UnitCell(a1,a2) 

b1 = addBasisSite!(uc, (0.0, 0.0))


M = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]


M4=ein"ij,kl->ijkl"(M,M)

addInteraction!(uc, [b1, b1, b1, b1], [(0, 0),(0, 0),(1, 0),(1, 0)],M4) 


ll=Lattice(uc, (2,2))
