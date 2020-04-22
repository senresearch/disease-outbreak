using LsqFit
using Optim
using DataFrames
import StatsBase.fitted

include("dynamics.jl")

struct SIRX <: Dynamics
    α::Float64
    β::Float64
    κ::Float64
    κ0::Float64
end

# default constructor
SIRX() = SIRX(1.0,1.0,0.0,0.0)

# another constructor
function getParams(κ::Float64,κ0::Float64,
    R0Free::Float64,TInfected::Float64,d::Dynamics)
    β = 1.0/TInfected
    α = R0Free*β
    return SIRX(α,β,κ,κ0)
end

# number of states
function nstates(d::SIRX)
    return 4
end

# names of states
function stateNames(d::SIRX)
    return ["S" "I" "R" "X"]
end

function R0Eff(d::SIRX)
    return d.α/(d.β+d.κ+d.κ0)
end

# change in a day
"""
change(s::Vector{Float64},d::SIRX)

Return the change in the state of the population in a day

s = state vector (S,I,R,X)
d = SIRX dynamics parameters
"""
function change(s::Vector{Float64},d::SIRX)
    # this is done purely for readability of the formula
    st = NamedTuple{(:S,:I,:R,:X)}(s)
    S = -( d.α * st.I + d.κ0) * st.S
    I = d.α*st.I*st.S - d.β*st.I - d.κ0*st.I - d.κ*st.I
    R = d.β*st.I + d.κ0*st.S
    X = (d.κ+d.κ0) * st.I
    return [S,I,R,X]
end

# initializer for states
function initialize(N::Float64,C0::Float64,IXRatio::Float64,d::SIRX)
    state = zeros(nstates(d))
    state[4] = C0/N # X0, page 1, Supplement
    state[2] = IXRatio*(C0/N) # I0, page 1, supplement
    state[1] = 1.0-state[2]-state[4] # S0
    return state
end

"""
caseModel(nt::Int64,N::Float64,C0::Float64,IXRatio::Float64,d::SIRX)

Returns vector of cases given model dynamics and population parameters.
"""
function caseModel(nt::Int64,N::Float64,C0::Float64,IXRatio::Float64,
    d::SIRX)
    s0 = initialize(N,C0,IXRatio,d)
    (d,ds) =evolve(N,s0,d,nt)
    return N.*d[4,:]
end

function estimatedStates(nt::Int64,N::Float64,C0::Float64,
    IXRatio::Float64,d::SIRX)
    s0 = initialize(N,C0,IXRatio,d)
    (d,ds) =evolve(N,s0,d,nt)
    s1 = DataFrame(N.*d',[:S,:I,:R,:X])
    return s1
end

function logCaseModelUnknown(nt::Int64,p::Vector{Float64},
    N::Float64,C0::Float64,R0Free::Float64,TInfected::Float64)
    d = getParams(p[1],p[2],R0Free,TInfected,SIRX())
    return log.(caseModel(nt,N,C0,p[3],d))
end

function sqrtCaseModelUnknown(nt::Int64,p::Vector{Float64},
    N::Float64,C0::Float64,R0Free::Float64,TInfected::Float64)
    d = getParams(p[1],p[2],R0Free,TInfected,SIRX())
    return sqrt.(caseModel(nt,N,C0,p[3],d))
end
