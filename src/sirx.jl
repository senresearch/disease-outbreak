include("dynamics.jl")

using LsqFit
using Optim

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

#######################

struct SIRXPopulation <: Population
    N::Float64
    C0::Float64
    IXRatio::Float64
end

# initializer for states
function initialize(pop::SIRXPopulation,d::SIRX)
    state = zeros(nstates(d))
    state[4] = pop.C0/pop.N # X0, page 1, Supplement
    state[2] = IXRatio*(pop.C0/pop.N) # I0, page 1, supplement
    state[1] = 1.0-state[2]-state[4] # S0
    return state
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

"""
caseModel(t::Vector{Float64},params::Vector{Float64},
    N::Float64,C0::Float64,R0::Float64,TI::Float64)

Returns vector of cases given model dynamics and population parameters.
"""
function caseModel(t::Vector{Float64},pop::SIRXPopulation,d::SIRX)
    s0 = initialize(pop,d)
    (d,ds) =evolve(pop.N,s0,d,length(t))
    return N.*d[4,:]
end

function logCaseModelUnknown(t::Vector{Float64},p::Vector{Float64},
    N::Float64,C0::Float64,R0Free::Float64,TInfected::Float64)
    d = getParams(p[1],p[2],R0Free,TInfected,SIRX())
    pop = SIRXPopulation(N,C0,p[3],d)
    return log.(caseModel(t,pop,d)
end

function rss(y::Vector{Float64},x::Vector{Float64},p::Vector{Float64})
    return sum((y-model(x,p)).^2)
end

x = convert(Vector{Float64},(1:length(hubei.ConfirmedCases)))
y = log.(hubei.ConfirmedCases)
fit = curve_fit(model,x,y,[log(0.1),log(0.1),log(1)])

fit
exp.(fit.param)
