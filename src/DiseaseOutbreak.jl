module DiseaseOutbreak

include("distancing.jl")
export makeCitizens, infectCitizens!, followCitizens!
export uninfected, newInfections, cumulativeInfections, infectionSpread

include("dynamics.jl")
export Dynamics, Population

include("evolution.jl")
export evolve

include("sirx.jl")
include("sei3r.jl")
export SEI3R
export SIRX
export getParams, initialize
export caseModel

include("casemodelfit.jl")
export CaseModelFitResult
export fitCaseModel, summary, fitted, estimatedStates

include("sirplot.jl")
export plotEvolution, plotFit
export pyplotEvolution, pyplotFit
end # module
