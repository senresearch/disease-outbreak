using Plots, Statistics, StatsBase
include("sirx.jl")
include("seir.jl")

Dynamics = Union{SEI3RDynamics,SIRXDynamics}

function evolve(state::Vector{Float64},d::Dynamics,ntime::Int64)
    states = zeros(nstates(d),ntime)
    deltaStates = zeros(nstates(d),ntime)
    states[:,1] = state
    for i in 2:ntime
        deltaStates[:,i] = change(states[:,i-1],d)
        states[:,i] = states[:,i-1].+deltaStates[:,i]
    end
    return (states=states,deltaStates=deltaStates)
end

function plotEvolution(states::Matrix{Float64},d::Dynamics)
    Plots.plot(Matrix(states'),label=stateNames(d))
end
