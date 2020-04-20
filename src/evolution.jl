using Plots, Statistics, StatsBase

include("sei3r.jl")
include("sirx.jl")

"""
s, ds = evolve(N::Float64,state::Vector{Float64},d::Dynamics,
    ntime::Int64)

Return a matrix of states over time (d) and corresponding changes in states per day
(ds) up to ntime.

N = the population size
state = a vector of initial state. See the function, 'initialze'.
d = a type of Dynamics
ntime = the length of time to evolve

"""

function evolve(N::Float64,state::Vector{Float64},d::Dynamics,
    ntime::Int64)
    states = zeros(nstates(d),ntime)
    deltaStates = zeros(nstates(d),ntime)
    states[:,1] = state
    for i in 2:ntime
        deltaStates[:,i] = change(states[:,i-1],d)
        states[:,i] = states[:,i-1].+deltaStates[:,i]
    end
    return (states=states,deltaStates=deltaStates)
end

"""
plotEvolution(N::Float64,states::Matrix{Float64},
    d::Dynamics;log10::Bool=false)

Generate a plot of states evolved for the SIR model.

N = the population size
states = a matrix of states evoloved over time. See the function, 'evolve'.
d = a type of Dynamics of the SIR model
log10 = a keyword argument (default is false) of scaling state values by log10

"""

function plotEvolution(N::Float64,states::Matrix{Float64},
    d::Dynamics;log10::Bool=false)
    if (log10==false)
        Plots.plot(N*Matrix(states'),label=stateNames(d))
    else
        Plots.plot(N*Matrix(states'),label=stateNames(d),ylab=:log10)
    end
end
