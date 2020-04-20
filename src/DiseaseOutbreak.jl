module DiseaseOutbreak

include("distancing.jl")
export makeCitizens, infectCitizens!, followCitizens!
export uninfected, newInfections, cumulativeInfections, infectionSpread

include("dynamics.jl")
export Dynamics, Population

# include("evolution.jl")
# export evolve, plotEvolution

include("sei3r.jl")
export SEI3R
export getParams, initialize

include("sirx.jl")
export SIRX
export getParams, initialize


end # module
