## implement Markov Chain with SEIR

abstract type Dynamics end
abstract type SEI3RDynamics <: Dynamics end

struct SEI3RDynamics
    α::Float64
    β::Vector{Float64}
    γ::Vector{Float64}
    p::Vector{Float64}
    μ::Float64
end


function vec2state(v::Vector{Float64})
    return State(S=v[1],E=v[2],I=copy(v[3:5]),R=v[6],D=v[7])
end

function oneDay(s::State,par::Dynamics)
    (S,E,I,R,D)=model0(s.S,s.E,s.I,s.R,s.D,
                       par.α,par.β,par.γ,par.p,par.μ)
    return State(S,E,I,R,D)
end

function change(s::Vector{Float64},d::SEI3RDynamics)
    ds = NamedTuple{(:S,:E,:I1,:I2,:I3,:R,:D)}(s)
    st = NamedTuple{(:S,:E,:I1,:I2,:I3,:R,:D)}(s)

    ds.S = -( d.β[1]*st.I[1] - d.β[2]*st.I[2] - d.β[3]*stI[3] ) * st.S
    ds.E =  ( s.β[1]*st.I[1] + d.β[2]*st.I[2] + s.β[3]*st.I[3] ) - d.α*st.E
    ds.I1 = d.α*st.E - (d.γ[1]+d.p[1]) * st.I[1]
    ds.I2 = d.p[1]*st.I[1] - (d.γ[2]+d.p[2]) * st.I[2]
    ds.I3 = d.p[2]*st.I[2] - (d.[3]+d.μ) * st.I[3]
    ds.R = d.γ[1]*st.I[1] + d.γ[2]*st.I[2] + d.γ[3]*st.I[3]
    ds.D = d.μ*st.I[3]
    return values(ds)
end


function update(s::State,d::Dynamics)
    return
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
