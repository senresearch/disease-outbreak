module DiseaseOutbreak

include("distancing.jl")
export makeCitizens, infectCitizens!, followCitizens!
export uninfected, newInfections, cumulativeInfections, infectionSpread

include("dynamics.jl")
export Dynamics, Population

include("evolution.jl")
export evolve, plotEvolution

include("sirx.jl")
include("sei3r.jl")
export SEI3R
export SIRX
export getParams, initialize
export fitCaseModel, caseModel, fitted



end # module
