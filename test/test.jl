using ClassicalSpinMonteCarlo
using LinearAlgebra

a1 = (1.0, )
uc = UnitCell(a1)

b1 = addBasisSite!(uc, (0.0,))
#b2 = addBasisSite!(uc, (1.0,))

M = -[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.5]
addInteraction!(uc, [b1, b1], [(0,),(-1,)],M)
#addInteraction!(uc, [b1, b2], [(0,),(0,)],M)
#addInteraction!(uc, [b1], [(0,)],[0.0,0.0,1])  


L = (16,)
lattice = Lattice(uc, L)

thermalizationSweeps = 50000
measurementSweeps = 10
beta = 10.0
m = MonteCarlo(lattice, beta, thermalizationSweeps, measurementSweeps)
run!(m)

sd=SpinDynamic(lattice,beta,20,nenergy=2000,denergy=0.05)


run!(sd)


sd1=SpinDynamic(sd.lattice,beta,10,nenergy=2000,denergy=0.05)

LLGrun!(sd1)

w,k,sp=Spectra(sd1)

ks=[[i] for i in k[1]]
w,sp1=Spectra(sd,ks)

k,s=structure_factor(sd,16)

k,s=structure_factor(m.lattice,16)

using Plots
heatmap(k[1],w,real.(sp[1,:,:]+sp[2,:,:])',yrange=(0.1,5),clims=(0,0.0001))#,aspect_ratio=1,xrange=(-2pi,2pi),yrange=(0,5),clims=(0,1),xlabel="kx", ylabel="ky")

heatmap(k[1],w,real.(sp1[1,:,:]+sp1[2,:,:]+sp1[3,:,:])',yrange=(0,5),clims=(0,0.0001))
heatmap(k[1],w,real.(sp1[3,:,:])',yrange=(0,5),clims=(0,0.001))