using LinearAlgebra
using FFTW

mutable struct SpinDynamic{T<:Lattice}
    lattice::T
	beta::Float64
	step::Float64
	num::Int
	time::Int
	samples::Vector{Array{Float64,3}}
	SpinDynamic(lattice::Lattice,beta::Float64,step::Float64,num::Int,time::Int) = new{typeof(lattice)}(lattice,beta,step,num,time,Array{Float64,3}[])
end

SpinDynamic(lattice::Lattice,beta::Float64,num::Int;nenergy::Int,denergy::Float64) = SpinDynamic(lattice,beta,2*pi/nenergy/denergy,num,nenergy)


function Spectra(sd::SpinDynamic{T}) where T<:Lattice
	L=sd.lattice.size
	k=map(x->0:2*pi/x:2*pi-2*pi/x,L)
	w=range(0,length=sd.time,step=2*pi/sd.time/sd.step)
	basis=sd.lattice.unitcell.basis
	primitive=sd.lattice.unitcell.primitive
	#spins=[reshape(s,3,length(basis),L...,sd.time) for s in sd.samples]
	ncolon=ntuple(x->Colon(),length(L)+1)
	spinkw=zeros(ComplexF64,3,length(basis),length(basis),L...,sd.time)
	for s in sd.samples
		sk=Array{ComplexF64,3+length(L)}(undef,3,length(basis),L...,sd.time)
		for i=1:length(basis),j=1:3
			sk[j,i,ncolon...]=fft(reshape(s,3,length(basis),L...,sd.time)[j,i,ncolon...])
		end
		for i=1:length(basis),j=i+1:length(basis),x=1:3
			spinkw[x,i,j,ncolon...]+=conj.(sk[x,i,ncolon...]).*sk[x,j,ncolon...]
		end
		for i=1:length(basis),x=1:3
			spinkw[x,i,i,ncolon...]+=abs.(sk[x,i,ncolon...]).^2
		end		
	end
	
	out=zeros(ComplexF64,3,L...,sd.time)
		for i=1:length(basis),j=i+1:length(basis)
			r=basis[i].-basis[j]
			r=[dot(r,p)./norm(p)^2 for p in primitive]
			ph=[exp(im*dot(ki,r)) for ki in Iterators.product(k...)]
		
			for x=1:3
				for y=1:sd.time
					spinkw[x,i,j,ncolon[1:length(L)]...,y]=spinkw[x,i,j,ncolon[1:length(L)]...,y].*ph
				end
				
				out[x,ncolon...]+=2*real.(spinkw[x,i,j,ncolon...])
			end
		end
		for i=1:length(basis),x=1:3
			out[x,ncolon...]+=spinkw[x,i,i,ncolon...]
		end	

	
	return w,k,real.(out)./(sd.num*(sd.lattice.length*sd.time)^2)

end

function Spectra(sd::SpinDynamic{T},ks::Vector{Vector{Float64}}) where T<:Lattice
	#L=sd.lattice.size
	w=range(0,length=sd.time,step=2*pi/sd.time/sd.step)
	basis=sd.lattice.unitcell.basis
	primitive=sd.lattice.unitcell.primitive
	positions=sd.lattice.sitePositions
	
	spinkw=zeros(ComplexF64,3,length(ks),sd.time)
	

	for (ik, k) in enumerate(ks)
		expqr=[exp(-im*dot(k,p)) for p in positions]
		for s in sd.samples
			for x=1:3
				Aq=Array{ComplexF64,1}(undef,sd.time)
			    for t=1:sd.time
					Aq[t]=sum(s[x,:,t].*expqr)
				end
				spinkw[x,ik,:]+=abs.(fft(Aq)).^2	
			end
			
		end
	end
	

	
	return w,real.(spinkw)./(sd.num*(sd.lattice.length*sd.time)^2)

end


function structure_factor(sd::SpinDynamic{T},Nk::Int) where T<:Lattice
	#L=sd.lattice.size

	basis=sd.lattice.unitcell.basis
	primitive=sd.lattice.unitcell.primitive
	positions=sd.lattice.sitePositions
	D=sd.lattice.size|>length
	
	nk=ntuple(i->Nk+1,D)
	spinkw=zeros(3,nk...)

    k = [collect(range(-2pi,2pi,length=Nk+1)) for _=1:D]
    #ky = collect(range(-2pi,2pi,length=Nk+1))
	
	for site in Iterators.product(range.(1,nk)...)

			ki=[k[i][site[i]] for i=1:D]
			expqr=[exp(-im*dot(ki,p)) for p in positions]
			for s in sd.samples
				for x=1:3

					spinkw[x,site...]+=abs(sum(s[x,:,1].*expqr))^2	
				end
				
			end

	end
	

	
	return k,spinkw./(sd.num*(sd.lattice.length)^2)

end


function run!(sd::SpinDynamic{T};method::Function=RK4) where T<:Lattice

    for i=1:sd.num
		for _ in 1:5*length(sd.lattice)
			metropolis(sd.lattice, sd.beta)
		end

        push!(sd.samples,cat(sd.lattice.spins,dims=3))
	end
	
	
	dynamic!(sd,method)

	
end

function LLGrun!(sd::SpinDynamic{T};method::Function=RK4) where T<:Lattice

    push!(sd.samples,cat(sd.lattice.spins,dims=3))
	for i=1:sd.num-1
		for _ in 1:length(sd.lattice)
			metropolis(sd.lattice, sd.beta)
		end

        push!(sd.samples,cat(sd.lattice.spins,dims=3))
	end
	
	
	dynamic!(sd,method)

	
end

function dynamic!(sd::SpinDynamic{T},method::Function=RK4) where T<:Lattice
	
	for i=1:sd.num
		 sd.samples[i]=method(sd.lattice.siteinteraction,sd.samples[i][:,:,1],sd.step,sd.time)
	end
	
end


function RK4ts(F::Matrix{Vector{Interaction}} ,first::Matrix{Float64},step::Float64,time::Int)
	step2=step/2
	step6=step/6
	out=[first]
	tsm=_ts(first)
	m0=first
	for _=2:time
		
		k1=_dfts(tsm,m0,F)
		k2=_dfts(tsm+step2*k1,F)
		k3=_dfts(tsm+step2*k2,F)
		k4=_dfts(tsm+step*k3,F)
			
		
		tsm+=step6*(k1+2*(k2+k3)+k4)
		m0=_inversts(tsm)
		push!(out,m0)
	end
	return cat(out...,dims=3)
end

function RK4(F::Matrix{Vector{Interaction}} ,first::Matrix{Float64},step::Float64,time::Int)
	step2=step/2
	step6=step/6
	out=[first]
	m=first
	for _=2:time
		
		k1=_df(m,F)
		k2=_df(m+step2*k1,F)
		k3=_df(m+step2*k2,F)
		k4=_df(m+step*k3,F)
			
		
		m+=step6*(k1+2*(k2+k3)+k4)
		
		for i=1:size(m,2)
			m[:,i]=normalize(m[:,i])
		end
		push!(out,m)
	end
	return cat(out...,dims=3)
end

function _df(y::Matrix{Float64},F::Matrix{Vector{Interaction}})
	ff=[sum(map(x->field(x,y),f)) for f in F]
	k=Matrix{Float64}(undef,size(y))
	for i=1:size(y,2)
		k[:,i]=cross(ff[:,i],y[:,i])
	end

	return k
end



function _dfts(y::Matrix{Float64},F::Matrix{Vector{Interaction}})
	y0=_inversts(y)
	ff=[sum(map(x->field(x,y0),f)) for f in F]
	k=Matrix{Float64}(undef,size(y))
	for i=1:size(y,2)
		s=tan(y[2,i])
		cosy=cos(y[1,i])
		siny=sin(y[1,i])
		k[1,i]=-(ff[1,i]*cosy+ff[2,i]*siny)/s+ff[3,i]
		k[2,i]=ff[2,i]*cosy-ff[1,i]*siny
		
	end
	return k
end


function _dfts(y::Matrix{Float64},y0::Matrix{Float64},F::Matrix{Vector{Interaction}})
	#y0=_inversts(y)
	ff=[sum(map(x->field(x,y0),f)) for f in F]
	k=Matrix{Float64}(undef,size(y))
	for i=1:size(y,2)
		s=tan(y[2,i])
		cosy=cos(y[1,i])
		siny=sin(y[1,i])
		k[1,i]=-(ff[1,i]*cosy+ff[2,i]*siny)/s+ff[3,i]
		k[2,i]=ff[2,i]*cosy-ff[1,i]*siny
		
	end
	return k
end


function _ts(spins::Matrix{Float64})
	out=Matrix{Float64}(undef,2,size(spins,2))
	out[2,:]=acos.(spins[3,:])
	out[1,:]=atan.(spins[2,:]./spins[1,:]).-(sign.(spins[1,:]).-1)*(pi/2)
	return out

end

function _inversts(spins::Matrix{Float64})
	out=Matrix{Float64}(undef,3,size(spins,2))
	out[3,:]=cos.(spins[2,:])
	sp=sin.(spins[2,:])
	out[1,:]=sp.*cos.(spins[1,:])
	out[2,:]=sp.*sin.(spins[1,:])
	return out

end



