using Random
using LinearAlgebra

function uniformOnSphere(rng = Random.GLOBAL_RNG)
    phi = 2.0 * pi * rand(rng)
    z = 2.0 * rand(rng) - 1.0;
    r = sqrt(1.0 - z * z)
    return (r * cos(phi), r * sin(phi), z)
end
function uniformOnSphere(site::Int, lattice::Lattice{D}, rng = Random.GLOBAL_RNG) where {D}
    s0 = getSpin(lattice, site)
    rspin = uniformOnSphere(rng)
    srs = (s0[1]*rspin[1] + s0[2]*rspin[2] + s0[3]*rspin[3] )
    return normalize([s0[1] - 2*rspin[1]*srs, s0[2] - 2*rspin[2]*srs, s0[3] - 2*rspin[3]*srs])  
end
function zeroheatbath_update0(site::Int, lattice::Lattice{D}) where {D}
    h = localH(lattice, site)
    s0 = getSpin(lattice, site)
    v = [s0[i] - 2*h[i] for i in 1:3]
    vn = norm(v)
    return v[1]/vn, v[2]/vn, v[3]/vn
end
function zeroheatbath_update(site::Int, lattice::Lattice{D}) where {D}
    h = localH(lattice, site)
    v = collect(Float64, h)
    res = -normalize(v)
    return res[1], res[2], res[3]
end
function overrelaxation_update(site::Int, lattice::Lattice{D}) where {D}
    h = localH(lattice, site)
    hn = norm(collect(Float64, h))
    s0 = getSpin(lattice, site)
    scale = 2*(s0[1]*h[1] + s0[2]*h[2] + s0[3]*h[3])/(hn^2)
    return normalize([scale*h[1]-s0[1], scale*h[2]-s0[2], scale*h[3]-s0[3]])
end
function heatbath_update(site::Int, lattice::Lattice{D}, beta::Float64, rng = Random.GLOBAL_RNG) where {D}
    h = localH(lattice, site)
    temp = collect(Float64, h)
    hn = norm(temp)#sqrt(h[1]^2 + h[2]^2 + h[3]^2)
    hv = normalize(temp)#[h[1]/hn, h[2]/hn, h[3]/hn]
    r = rand(rng)
    ϕ = 2.0 * pi * rand(rng)
    bg = (hn*beta)
    if abs(bg) > 200.0 
        z = -sign(bg)
    else
        z = -(log(exp(-bg)*(1-r) + r*exp(bg)))/(bg)  
    end
    θ = acos(z)
    newspin_h = Float64[sin(θ)*cos(ϕ) sin(θ)*sin(ϕ) z]
    RT = _euleranglesT(hv)
    newspin = newspin_h*RT
    return newspin[1], newspin[2], newspin[3]
end

function _euleranglesT( v2::Vector{Float64} ) ::Matrix{Float64}
    # v1 = [0 0 1]
    v₁ = Float64[0, 0, 1.0]
    v₂ = v2
    u = Float64[-v₂[2], v₂[1], 0]#cross(v₁, v₂)
    M = Float64[v₁ v₂ u]
    MT = transpose(M)#Float64[0 0 1; v₂[1] v₂[2] v₂[3]; u[1] u[2] u[3]]
    vv = v₂[3]
    # V = [0 -1 0; 1 0 0; 0 0 vv==-1 ? 0 : 1/(1+vv)]
    Vt = [0 1 0; -1 0 0; 0 0 vv==-1 ? 0 : 1/(1+vv)]
    RT = vv*I(3) + M*Vt*MT
    return RT
end



#=
function _localH(M::InteractionMatrix, s2)::NTuple{3, Float64}
    return  (M.m11 * s2[1] + M.m12 * s2[2] + M.m13 * s2[3]),  (M.m21 * s2[1] + M.m22 * s2[2] + M.m23 * s2[3]) , (M.m31 * s2[1] + M.m32 * s2[2] + M.m33 * s2[3])   
end
=#
"""
    localH(lattice::Lattice{D}, site::Int)::NTuple{3, Float64}  where {D}

only for no on-site interaction. e.g., J*s_i*s_i
"""
function field(i::Interaction{R},spins::Matrix{Float64}) where R
	return R*i.value*prod(spins[i.indices])
end



function localH(lattice::Lattice{D}, site::Int)::NTuple{3, Float64}  where {D}
    #two-spin interactions
    x,y,z = lattice.siteinteraction[:,site]

    hx, hy, hz = 0.0, 0.0, 0.0
    for i in x 
        hx += field(i,lattice.spins)
    end
	for i in y 
        hy += field(i,lattice.spins)
    end
	for i in z 
        hz += field(i,lattice.spins)
    end

    return hx , hy , hz 
end


#=
function exchangeEnergy(s1, M::InteractionMatrix, s2)::Float64
    return s1[1] * (M.m11 * s2[1] + M.m12 * s2[2] + M.m13 * s2[3]) + s1[2] * (M.m21 * s2[1] + M.m22 * s2[2] + M.m23 * s2[3]) + s1[3] * (M.m31 * s2[1] + M.m32 * s2[2] + M.m33 * s2[3])
end
=#

function getEnergy(lattice::Lattice{D})::Float64 where {D}
    energy = 0.0

    for (key,value) in pairs(lattice.interactions)
		
        energy += value*prod(lattice.spins[key])
    end

    return energy
end

function energy(i::Interaction{R},spins::Matrix{Float64},newState::Float64) where R
	ds=newState^R-spins[i.index]^R
	return ds*i.value*prod(spins[i.indices])
end

function getEnergyDifference(lattice::Lattice{D}, site::Int, newState::Tuple{Float64,Float64,Float64})::Float64 where {D}
    dE = 0.0
	
    x,y,z = lattice.siteinteraction[:,site]

    for i in x 
        dE += energy(i,lattice.spins,newState[1])
    end
	for i in y 
        dE += energy(i,lattice.spins,newState[2])
    end
	for i in z 
        dE += energy(i,lattice.spins,newState[3])
    end
    return dE
end

function getMagnetization(lattice::Lattice{D}) where {D}

	#=
    mx, my, mz = 0.0, 0.0, 0.0
    for i in 1:length(lattice)
        spin = getSpin(lattice, i)
        mx += spin[1]
        my += spin[2]
        mz += spin[3]
    end
	=#
    return sum(lattice.spins,dims=2)[:] / length(lattice)
end

function getCorrelation(lattice::Lattice{D}, basis::Bool=true) where {D}
    nsite = basis == true ? length(lattice.unitcell.basis) : length(lattice)
    corr = zeros(length(lattice), nsite)
    if basis == false
        for i in 1:nsite
            s0 = getSpin(lattice, i)
            for j in i:length(lattice)
                corr[j,i] = dot(s0, getSpin(lattice, j))
            end
        end
        return Symmetric(corr, :L)
    else
        for i in 1:nsite
            s0 = getSpin(lattice, i)
            for j in 1:length(lattice)
                corr[j,i] = dot(s0, getSpin(lattice, j))
            end
        end
        return corr
    end 
end
# function getCorrelation(lattice::Lattice{D,N}) where {D,N}
#     corr = zeros(length(lattice), length(lattice))
#     for i in 1:length(lattice)
#         s0 = getSpin(lattice, i)
#         for j in 1:length(lattice)
#             corr[j,i] = dot(s0, getSpin(lattice, j))
#         end
#     end
#     return corr
# end