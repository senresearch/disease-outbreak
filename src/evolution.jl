using Plots, Statistics, StatsBase
include("sei3r.jl")
include("sirx.jl")

Dynamics = Union{SEI3RDynamics,SIRXDynamics}

function evolve(N::Float64,state::Vector{Float64},d::Dynamics,ntime::Int64)
    states = zeros(nstates(d),ntime)
    deltaStates = zeros(nstates(d),ntime)
    states[:,1] = state
    for i in 2:ntime
        deltaStates[:,i] = change(states[:,i-1],d)
        states[:,i] = states[:,i-1].+deltaStates[:,i]
    end
    return (states=states,deltaStates=deltaStates)
end

function plotEvolution(N::Float64,states::Matrix{Float64},d::Dynamics)
    Plots.plot(N*Matrix(states'),label=stateNames(d))
end
