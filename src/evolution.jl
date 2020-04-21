using Plots, Statistics, StatsBase


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
