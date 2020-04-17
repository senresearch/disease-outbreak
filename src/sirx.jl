include("dynamics.jl")

struct SIRX <: Dynamics
    α::Float64
    β::Float64
    κ::Float64
    κ0::Float64
end

# default constructor
SIRX() = SIRX(1.0,1.0,0.0,0.0)

"""
state0 = initialize(N::Float64,C0::Float64,IXRatio::Float64,
    d::SIRX)

Return a vector of inital state values.
"""

function initialize(N::Float64,C0::Float64,IXRatio::Float64,
    d::SIRX)
    state = zeros(nstates(d))
    state[4] = C0/N # X0, page 1, Supplement
    state[2] = IXRatio*(C0/N) # I0, page 1, supplement
    state[1] = 1.0-state[2]-state[4] # S0
    return state
end

# change in a day
"""
change(s::Vector{Float64},d::SIRX)

Retutn the change in the state of the population in a day.

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

function nstates(d::SIRX)
    return 4
end
function stateNames(d::SIRX)
    return ["S" "I" "R" "X"]
end

"""
param0 = getParams(κ::Float64,κ0::Float64,
    R0Free::Float64,TInfected::Float64,d::Dynamics)

Return inital parameter values of SIRX dynamics parameters.

κ = the rate of quarantine measures for symptomatic infected
    indviduals in (I)nfected state
κ0 = the containment rate effective in both (S)useptible and I states
      i.e. social distancing, curfews, etc.
R0Free = basic (unconstrained) reproduction number computed by α/β.
TInfected = the average time for an individual to remain infectious
     in I before (R)emoved
d = a type of Dynamics, SIRX

"""

function getParams(κ::Float64,κ0::Float64,
    R0Free::Float64,TInfected::Float64,d::Dynamics)
    β = 1.0/TInfected
    α = R0Free*β
    return SIRX(α,β,κ,κ0)
end
