
using RecipesBase: RecipesBase, @recipe, @series, @layout

function structure_factor(Nk::Int, lattice::Lattice{D}, correlation::Matrix{Float64}) where {D}
    kx = collect(range(-2pi, 2pi, length=Nk))
    ky = collect(range(-2pi, 2pi, length=Nk))
    structurefactor = zeros(ComplexF64,Nk,Nk)
    nsite, n2 = size(correlation, 2), size(correlation,1)
    for i in 1:Nk
        for j in 1:Nk
            z = 0.0
            # Compute Fourier transformation at momentum (kx, ky). The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
            for b in 1:nsite
                for k in 1:n2
                    z += cos(dot((kx[i],ky[j]), getSitePosition(lattice,k).-getSitePosition(lattice,b))) * correlation[k,b]
                end
            end
            structurefactor[j,i] = z / (length(lattice) * nsite)
        end
    end
        return kx, ky, structurefactor
end

function structure_factor(Nk::Int, lattice::Lattice{D}) where {D}
    kx = collect(range(-2pi,2pi,length=Nk))
    ky = collect(range(-2pi,2pi,length=Nk))
    structurefactor = zeros(ComplexF64,Nk,Nk)
    for i in 1:Nk
        for j in 1:Nk
            z = 0.0
        # Compute Fourier transformation at momentum (kx, ky). The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
            for b in 1:length(lattice)
                for k in 1:length(lattice)
                    correlation = dot(lattice.spins[:,k], lattice.spins[:,b])
                    z += exp(1im*dot((kx[i],ky[j]),getSitePosition(lattice,k).-getSitePosition(lattice,b))) * correlation
                end
            end
            structurefactor[j,i] = z / (length(lattice) )
        end
    end
    return kx,ky,structurefactor
end

function structure_factor(lattice::Lattice{D},Nk::Int) where {D}

    nk=ntuple(i->Nk+1,D)
	k = [collect(range(-2pi,2pi,length=Nk+1)) for _=1:D]
	#kx = collect(range(-2pi,2pi,length=Nk+1))
    #ky = collect(range(-2pi,2pi,length=Nk+1))
    structurefactor = zeros(3,nk...)
	positions=lattice.sitePositions
	spins=lattice.spins
    for site in Iterators.product(range.(1,nk)...)

			ki=[k[i][site[i]] for i=1:D]
			expqr=[exp(-im*dot(ki,p)) for p in positions]
			for x=1:3

				structurefactor[x,site...] = abs(sum(spins[x,:].*expqr))^2
			
			end			

    end
    return k,structurefactor./length(lattice)^2 
end

"""
    @recipe plot(lattice::Lattice{D}; spinscale=0.5) where {D}

Define the recipe for the visualization of a lattice and the spin configuration.
"""
@recipe function plot(lattice::Lattice{D}; spinscale=0.5,x::Int=1,y::Int=2) where {D}
    positions = zeros(2, length(lattice)) #NTuple{length(lattice), Float64}
    for i in 1:length(lattice)
        positions[:, i] = collect(getSitePosition(lattice, i))
    end
    spin_vector = zeros(3, length(lattice))
    for i in 1:length(lattice)
        spin_vector[:, i] = (normalize(lattice.spins[:,i])*spinscale)
    end
    titlefontsize --> 10
    legend := false
    aspect_ratio := :equal
    @series begin
        seriestype := :scatter
        markercolor --> :black
        positions[1, :], positions[2, :]
    end
    @series begin
        arrow --> true
        lw --> 1.5
        xp = [[positions[1, i] - spin_vector[x, i]/2, positions[1, i] + spin_vector[x, i]/2] for i in 1:lattice.length ]
        yp = [[positions[2, i] - spin_vector[y, i]/2, positions[2, i] + spin_vector[y, i]/2] for i in 1:lattice.length ]
        xp, yp
    end
end

# plt1=plot( positions[1,:],[positions[2,:]],legend=false,st=:scatter, grid=false, size=(400, 400))
# scale = 0.5
# spin_vector = zeros(3,length(lattice))
# for i in 1:length(lattice)
#     spin_vector[:, i] = normalize(m.lattice.spins[:,i])*scale
# end
# # Plot an arrow representing the spin vector
# quiver!(plt1,positions[1,:].-spin_vector[1,:]/2, positions[2,:].-spin_vector[2,:]/2, quiver=(spin_vector[1,:], spin_vector[2,:]), color=:blue, lw=2, label="Spin Vector")
# display(plt1)
# plt2=plot( positions[1,:],[positions[2,:]],legend=false,st=:scatter, grid=false, size=(400, 400))
# quiver!(plt2,positions[1,:].-spin_vector[1,:]/2, positions[2,:].-spin_vector[3,:]/2, quiver=(spin_vector[1,:], spin_vector[3,:]), color=:blue, lw=2, label="Spin Vector")
# display(plt2)