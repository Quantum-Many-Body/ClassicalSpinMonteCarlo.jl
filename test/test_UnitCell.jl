using OMEinsum
using ClassicalSpinMonteCarlo
using LinearAlgebra


a1 = (3/2, sqrt(3)/2)
a2 = (3/2, -sqrt(3)/2)
uc = UnitCell(a1,a2) 

b1 = addBasisSite!(uc, (0.0, 0.0))
b2 = addBasisSite!(uc, (1.0, 0.0)) 

M = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
addInteraction!(uc, [b1, b2], [(0, 0),(0, 0)],M) 
addInteraction!(uc, [b1, b2], [(0, 0),(-1, 0)],M) 
addInteraction!(uc, [b1, b2], [(0, 0),(0, -1)],M) 

addInteraction!(uc, [b1], [(0, 0)],[1.0,0.0,0.0]) 

M4=ein"ij,kl->ijkl"(M,M)

addInteraction!(uc, [b1, b1, b2, b2], [(0, 0),(0, 0),(0, 0),(0, 0)],M4) 
