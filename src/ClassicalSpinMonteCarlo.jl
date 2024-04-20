module ClassicalSpinMonteCarlo


include("UnitCell.jl")
export UnitCell, addInteraction!, addBasisSite! #,setInteractionOnsite!, setField!
#include("InteractionMatrix.jl")
include("Lattice.jl")
export Lattice, size, length, getSpin, setSpin!, getSitePosition

include("Observables.jl")
export Observables
include("Spin.jl")
export getEnergy, getMagnetization, getCorrelation

include("MonteCarlo.jl")
export MonteCarlo, run!

include("Helper.jl")
include("IO.jl")
export writeMonteCarlo, readMonteCarlo

include("Structurefactor.jl")
export structure_factor

include("SpinDynamic.jl")
export Spectra, SpinDynamic, run!,LLGrun!

using Reexport
@reexport using BinningAnalysis

end
