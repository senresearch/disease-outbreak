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

# initializer for states
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
function caseModel(t::Vector{Float64},params::Vector{Float64},
    N::Float64,C0::Float64,R0::Float64,TI::Float64)
    d = getParams(params[1],params[2],R0,TI,SIRX())
    s0 = initialize(N,C0,params[3],d)
    (d,ds) =evolve(N,s0,d,length(t))
    return N.*d[4,:]
end

function model(x::Vector{Float64},p::Vector{Float64})
    return log.(caseModel(x,exp.(p),hubeiPop,
                hubei.ConfirmedCases[1]*1.0,6.2,8.0))
end

function ss(y::Vector{Float64},x::Vector{Float64},p::Vector{Float64})
    return sum((y-model(x,p)).^2)
end

x = convert(Vector{Float64},(1:length(hubei.ConfirmedCases)))
y = log.(hubei.ConfirmedCases)
fit = curve_fit(model,x,y,[log(0.1),log(0.1),log(1)])

fit
exp.(fit.param)
