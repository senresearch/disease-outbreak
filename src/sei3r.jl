## implement Markov Chain with SEIR
include("dynamics.jl")

struct SEI3R <: Dynamics
    α::Float64
    β::Vector{Float64}
    γ::Vector{Float64}
    p::Vector{Float64}
    μ::Float64
end

# default constructor
SEI3R() = SEI3R(1.0,[1.0, 0.0, 0.0],zeros(3),zeros(2),0.0)

# change in a day
"""
change(s::Vector{Float64},d::SEI3R)

Retutn the change in the state of the population in a day

s = state vector (S,E,I1,I2,I3,R,D)
d = SEI3R dynamics parameters
"""
function change(s::Vector{Float64},d::SEI3R)
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

function initialize(E::Float64,d::SEI3R)
    state = zeros(nstates(d))
    state[2] = E
    state[1] = 1-E
    return state
end

function initialize(E::Float64,I1::Float64,d::SEI3R)
    state = zeros(nstates(d))
    state[3] = I1
    state[2] = E
    state[1] = 1-E-I1
    return state
end

function nstates(d::SEI3R)
    return 7
end
function stateNames(d::SEI3R)
    return ["S" "E" "I1" "I2" "I3" "R" "D"]
end


function getParams(IncubPeriod::Float64,DurMildInf::Float64,
    MildRate::Float64,SevereRate::Float64,CriticalRate::Float64,
    FracSevere::Float64,FracCritical::Float64,
    DurHosp::Float64,TimeICUStay::Float64,ICUDeathRate::Float64,
    d::SEI3R)

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

    return SEI3R(α,β,γ,p,μ)
end

function estimatedStates(nt::Int64,N::Int64,C0::Int64,
    CIRatio::Float64,d::SEI3R)
    s0 = initialize(0.0,(C0/N)/CIRatio,d)
    (d,ds) =evolve(N*1.0,s0,d,nt)
    s1 = DataFrame(N.*d',[:S,:E,:I1,:I2,:I3,:R,:D])
    return s1
end
