```@meta
CurrentModule = ClassicalSpinMonteCarlo
```
# Honeycomb Antiferromagnet

Spin dynamics of XXZ model in honeycomb lattice.

```@example
using ClassicalSpinMonteCarlo
using LinearAlgebra

#AB子格构造元胞
a1 = (1.0, 0.0)
a2=(1.0/2, √3/2)
uc = UnitCell(a1,a2)

#平移基矢
b1 = addBasisSite!(uc, (0.0,0.0))
b2 = addBasisSite!(uc, (0.5, √3/6))

#自旋相互作用矩阵(sx sy sz)*J*(sx sy sz)'
J1y = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.10]
J1x = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.10]
J1z = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.10]

#dxy= -0.23*1.08/3*[1.0 1.0 1.0;1.0 1.0 1.0;1.0 1.0 1.0]

#h=[1.0,1.0,1.0]

#JB= zeros(3,3,3,3)
#for i=1:3,j=1:3
#	JB[i,j,i,j]=-0.35
#end


#onsite相互作用
#addInteraction!(uc, [b1, b1], [(0,0),(0,0)],dxy)
#addInteraction!(uc, [b2, b2], [(0,0),(0,0)],dxy)

#addInteraction!(元胞, 子格[向量], 对应的正格矢[向量],相互作用张量[Array{n}]n=自旋算符数量)

addInteraction!(uc, [b1, b2], [(0,0),(0,0)],J1y)
addInteraction!(uc, [b1, b2], [(0,0),(-1,0)],J1x)
addInteraction!(uc, [b1, b2], [(0,0),(0,-1)],J1z)


#多极子相互作用
#addInteraction!(uc, [b1, b1,b2,b2], [(0,0),(0,0),(0,0),(0,0)],JB)
#addInteraction!(uc, [b1, b1,b2,b2], [(0,0),(0,0),(-1,0),(-1,0)],JB)
#addInteraction!(uc, [b1, b1,b2,b2], [(0,0),(0,0),(0,-1),(0,-1)],JB)


#外场
#addInteraction!(uc, [b1], [(0,)],[0.0,0.0,1])  

#退火轮数
thermalizationSweeps = 5000
#样本数目
measurementSweeps = 3
#1/温度
beta = 3.0#8.86817

#格子尺寸
L = (10,10)

#初始化格子	
lattice = Lattice(uc, L)
#可以初始化自旋构型Matrix(Float,3(自旋矢量xyz三个分量),格点数目)
#lattice.spins=spins


#初始化MC
m = MonteCarlo(lattice, beta, thermalizationSweeps, measurementSweeps)
#m = MonteCarlo(格子, 1/温度, 退火轮次, 样本数目)

#运行
run!(m)

#计算比热与误差
#c(e) = beta * beta * (e[2] - e[1] * e[1]) * length(m.lattice)
#∇c(e) = [-2.0 * beta * beta * e[1] * length(m.lattice), beta * beta * length(m.lattice)]
#heat[i] = mean(m.observables.energy, c)
#dheat[i] = std_error(m.observables.energy, ∇c)


#计算磁化与误差
#mean(m.observables.magnetization)
#dmagnetization[i] = std_error(m.observables.magnetization)

#计算自旋结构因子
k,s=structure_factor(m.lattice,30)
#倒空间[vector[kx],vector[ky]],自旋结构因子=structure_factor(m.lattice,30)

#heatmap(k[1],k[2],(s[1,:,:]+s[2,:,:]+s[3,:,:])')


   
#初始化自旋激发谱
sd=SpinDynamic(m.lattice,beta,3,nenergy=2000,denergy=0.02)
#sd=SpinDynamic(格子,1/温度,样本数目,能量数目,能量间隔)
#可选：sd=SpinDynamic(格子,1/温度,步长=2*pi/能量数目/能量间隔，样本数目,能量间隔)

#先退火5轮再计算
#run!(sd)


#初始化自旋
spin=zeros(3,2*prod(L))
spin[3,1:2:end].=1
spin[3,2:2:end].=-1
sd.lattice.spins=spin
#直接计算
LLGrun!(sd)

#计算结构因子
#k,s=structure_factor(sd,30)

#k空间高对称线
N=L[1]
ks=[[0,4*pi/sqrt(3)/N*i] for i=0:N]
ks1=[[2*pi/N*i,-2*pi/sqrt(3)/N*i] for i=0:N]
ks2=[[-2*pi/N*i,-2*pi/sqrt(3)/N*i] for i=0:N]
kp=[[2*pi/N*i,4*pi/sqrt(3)/N*(N-i)-2*pi/sqrt(3)/N*i] for i=1:N-1]
kp1=[[2*pi/N*(N-2*i),-2*pi/sqrt(3)] for i=1:N-1]
kp2=[[-2*pi/N*i,4*pi/sqrt(3)/N*(N-i)-2*pi/sqrt(3)/N*i] for i=1:N-1]
k=cat(ks,kp,ks1[end:-1:1],dims=1)
k1=cat(ks1,kp1,ks2[end:-1:1],dims=1)
k2=cat(ks,kp2,ks2[end:-1:1],dims=1)



#计算激发谱
w,sp=Spectra(sd,k)
#能量[向量],谱[Array(3[xx,yy,xx],k点数目,能量数目)]=Spectra(sd,倒空间点)
#可选：用FFT计算全空间激发谱：w,k,sp=Spectra(sd)





sp=log.(sp.+1)


using Plots
heatmap(1:3*N+1,w,(sp[1,:,:]+sp[2,:,:]+sp[3,:,:])',yrange=(0,4),clims=(0,maximum(sp)/10000))
```


