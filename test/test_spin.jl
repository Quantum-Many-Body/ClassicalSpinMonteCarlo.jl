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

#M4=ein"ij,kl->ijkl"(M,M)

#addInteraction!(uc, [b1, b1, b2, b2], [(0, 0),(0, 0),(0, 0),(0, 0)],M4) 


ll=Lattice(uc, (2,2))

a=getEnergy(ll)

a-(sum(ll.spins[:,1].*ll.spins[:,2]+ll.spins[:,3].*ll.spins[:,2]+ll.spins[:,5].*ll.spins[:,2] +ll.spins[:,3].*ll.spins[:,4]+ll.spins[:,7].*ll.spins[:,4]+ll.spins[:,1].*ll.spins[:,4] +ll.spins[:,5].*ll.spins[:,6]+ll.spins[:,7].*ll.spins[:,6]+ll.spins[:,1].*ll.spins[:,6] +ll.spins[:,7].*ll.spins[:,8]+ll.spins[:,3].*ll.spins[:,8]+ll.spins[:,5].*ll.spins[:,8]) +ll.spins[1,1]+ll.spins[1,3]+ll.spins[1,5]+ll.spins[1,7])

getMagnetization(ll)

getCorrelation(ll)


ClassicalSpinMonteCarlo.localH(ll,1).-(ll.spins[1,2]+ll.spins[1,4]+ll.spins[1,6]+1,ll.spins[2,2]+ll.spins[2,4]+ll.spins[2,6],ll.spins[3,2]+ll.spins[3,4]+ll.spins[3,6])

ns=rand(3)|>Tuple
ClassicalSpinMonteCarlo.getEnergyDifference(ll, 1,ns)-(ClassicalSpinMonteCarlo.localH(ll,1).*(ns.-ll.spins[:,1])|>sum)