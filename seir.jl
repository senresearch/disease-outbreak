## implement Markov Chain with SEIR

mutable struct State
    S::Float64
    E::Float64
    I::Vector{Float64}
    R::Float64
    D::Float64
end

struct Dynamics
    α::Float64
    β::Vector{Float64}
    p::Vector{Float64}
    γ::Vector{Float64}
    μ::Float64
end

function add(x::State,y::State)
    return State(x.S+y.S,x.E+y.E,x.I.+y.I,x.R+y.R,x.D+y.D)
end
function state2vec(s::State)
    return [ s.S s.E s.I s.R s.D ][1,:]
end

function vec2state(v::Vector{Float64})
    return State(S=v[1],E=v[2],I=copy(v[3:5]),R=v[6],D=v[7])
end

function oneDay(s::State,par::Dynamics)
    (S,E,I,R,D)=model0(s.S,s.E,s.I,s.R,s.D,
                       par.α,par.β,par.γ,par.p,par.μ)
    return State(S,E,I,R,D)
end

function model0(S,E,I,R,D,α,β,γ,p,μ)
    dS = -( β[1]*I[1] - β[2]*I[2] - β[3]*I[3] ) * S
    dE =  ( β[1]*I[1] + β[2]*I[2] + β[3]*I[3] ) - α*E
    dI=zeros(3)
    dI[1] = α*E - (γ[1]+p[1]) * I[1]
    dI[2] = p[1]*I[1] - (γ[2]+p[2]) * I[2]
    dI[3] = p[2]*I[2] - (γ[3]+μ) * I[3]
    dR = γ[1]*I[1] + γ[2]*I[2] + γ[3]*I[3]
    dD = μ*I[3]
    return S+dS, E+dE, I.+dI, R+dR, D+dD
end





function getParams(IncubPeriod::Float64,DurMildInf::Float64,
    MildRate::Float64,SevereRate::Float64,CriticalRate::Float64,
    FracMild::Float64,FracCritical::Float64,FracSevere::Float64,
    DurHosp::Float64,TimeICUDeath::Float64,CFR::Float64)

    β=[MildRate,SevereRate,CriticalRate]
    γ=zeros(3)
    p=zeros(2)


    α=1/IncubPeriod
    γ[1]=(1.0/DurMildInf)*FracMild
    p[1]=(1.0/DurMildInf)-γ[1]

    p[2]=(1.0/DurHosp)*(FracCritical/(FracSevere+FracCritical))
    γ[2]=(1.0/DurHosp)-p[2]

    μ=(1.0/TimeICUDeath)*(CFR/FracCritical)
    γ[3]=(1.0/TimeICUDeath)-μ


    return Dynamics(α,β,p,γ,μ)
end

d=getParams(5.0,6.0,0.5,0.1,0.1,0.15,0.05,0.1,6.0,8.0,2.0)

n = State(1000.0,0.0,zeros(3),0.0,0.0)

n = oneDay(n,d)
for i in 1:300
  n = oneDay(n,d)
end
