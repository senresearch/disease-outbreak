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

# function rss(y::Vector{Float64},x::Vector{Float64},p::Vector{Float64})
#    return sum((y-model(x,p)).^2)
# end

struct CaseModelFitResult
    κ::Float64
    κ0::Float64
    IXRatio::Float64
    inputs::NamedTuple
    fit::LsqFit.LsqFitResult
end

function fitCaseModel(nt::Int64,cases::Vector{Float64},
    N::Float64,R0Free::Float64,TInfected::Float64,
    p0::Vector{Float64})
    inputs = (nt=nt,C0=cases[1],N=N,R0Free=R0Free,TInfected=TInfected,C=cases)
    model(t,p) = logCaseModelUnknown(nt,exp.(p),N,cases[1],R0Free,TInfected)
    fit = curve_fit(model,(1:nt)*1.0,log.(cases),log.(p0))
    return CaseModelFitResult( exp(fit.param[1]), exp(fit.param[2]),
    exp(fit.param[3]), inputs, fit )
end


function summary(fit::CaseModelFitResult)
    r0eff = R0Eff(getParams(fit.κ,fit.κ0,fit.inputs.R0Free,
            fit.inputs.TInfected,SIRX()))
    res = return (κ=fit.κ,κ0=fit.κ0,IXRatio=fit.IXRatio,R0Eff=r0eff)
end

function fitted(fit::CaseModelFitResult,nt::Int64=0)
    if (nt==0)
        nt = fit.inputs.nt
    end
    return exp.(logCaseModelUnknown(nt,
                exp.(fit.fit.param), fit.inputs.N, fit.inputs.C0,
                fit.inputs.R0Free,fit.inputs.TInfected))
end

function estimatedStates(fit::CaseModelFitResult,nt::Int64=0)
    if (nt==0)
        nt = fit.inputs.nt
    end
    d = getParams(exp(fit.fit.param[1]),
                  exp(fit.fit.param[2]), fit.inputs.R0Free,
                  fit.inputs.TInfected, SIRX() )
    return estimatedStates( nt, fit.inputs.N, fit.inputs.C0,
           exp(fit.fit.param[3]), d )
end


function plotfit(fit::CaseModelFitResult,nt::Int64=0)
    if (nt==0)
        nt = fit.inputs.nt
    end
    plot(fit.inputs.C,yaxis=:log,seriestype=:scatter,
                color=:black,label="actual")
    plot!(estimatedStates(fit,nt)[:I],
                yaxis=:log,label="infected")
    plot!(fitted(fit,nt),yaxis=:log, label="fitted",color=:blue)
end
