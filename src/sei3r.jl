## implement Markov Chain with SEIR
include("dynamics.jl")

struct SEI3RDynamics <: Dynamics
    α::Float64
    β::Vector{Float64}
    γ::Vector{Float64}
    p::Vector{Float64}
    μ::Float64
end

# change in a day
"""
change(s::Vector{Float64},d::SEI3RDynamics)

Retutn the change in the state of the population in a day

s = state vector (S,E,I1,I2,I3,R,D)
d = SEI3R dynamics parameters
"""
function change(s::Vector{Float64},d::SEI3RDynamics)
    # this is done purely for readability of the formula
    st = NamedTuple{(:S,:E,:I1,:I2,:I3,:R,:D)}(s)
    # N = sum(s) # population size
    S = -( d.β[1]*st.I1 + d.β[2]*st.I2 + d.β[3]*st.I3 )*st.S
    E =  ( d.β[1]*st.I1 + d.β[2]*st.I2 + d.β[3]*st.I3 )*st.S - d.α*st.E
    I1 = d.α*st.E - (d.γ[1]+d.p[1]) * st.I1
    I2 = d.p[1]*st.I1 - (d.γ[2]+d.p[2]) * st.I2
    I3 = d.p[2]*st.I2 - (d.γ[3]+d.μ) * st.I3
    R = d.γ[1]*st.I1 + d.γ[2]*st.I2 + d.γ[3]*st.I3
    D = d.μ*st.I3
    return [S,E,I1,I2,I3,R,D]
end

function nstates(d::SEI3RDynamics)
    return 7
end
function stateNames(d::SEI3RDynamics)
    return ["S" "E" "I1" "I2" "I3" "R" "D"]
end

function getParams(IncubPeriod::Float64,DurMildInf::Float64,
    MildRate::Float64,SevereRate::Float64,CriticalRate::Float64,
    FracSevere::Float64,FracCritical::Float64,
    DurHosp::Float64,TimeICUStay::Float64,ICUDeathRate::Float64)

    CFR = ICUDeathRate*FracCritical
    FracMild = 1.0-FracSevere-FracCritical

    β=[MildRate,SevereRate,CriticalRate]
    γ=zeros(3)
    p=zeros(2)

    α=1/IncubPeriod
    γ[1]=(1.0/DurMildInf)*FracMild
    p[1]=(1.0/DurMildInf)-γ[1]

    p[2]=(1.0/DurHosp)*(FracCritical/(FracSevere+FracCritical))
    γ[2]=(1.0/DurHosp)-p[2]

    μ=(1.0/TimeICUStay)*(CFR/FracCritical)
    γ[3]=(1.0/TimeICUStay)-μ

    return SEI3RDynamics(α,β,γ,p,μ)
end

function initialize(E::Float64,d::SEI3RDynamics)
    state = zeros(nstates(d))
    state[2] = E
    state[1] = 1-E
    return state
end
