## implement Markov Chain with SEIR

# abstract type Dynamics end
# abstract type SEI3RDynamics <: Dynamics end

struct SEI3RDynamics
    α::Float64
    β::Vector{Float64}
    γ::Vector{Float64}
    p::Vector{Float64}
    μ::Float64
end


function change(s::Vector{Float64},d::SEI3RDynamics)
    # this is done purely for readability of the formula
    st = NamedTuple{(:S,:E,:I1,:I2,:I3,:R,:D)}(s)

    S = -( d.β[1]*st.I1 - d.β[2]*st.I2 - d.β[3]*st.I3 ) * st.S
    E =  ( d.β[1]*st.I1 + d.β[2]*st.I2 + d.β[3]*st.I3 ) - d.α*st.E
    I1 = d.α*st.E - (d.γ[1]+d.p[1]) * st.I1
    I2 = d.p[1]*st.I1 - (d.γ[2]+d.p[2]) * st.I2
    I3 = d.p[2]*st.I2 - (d.γ[3]+d.μ) * st.I3
    R = d.γ[1]*st.I1 + d.γ[2]*st.I2 + d.γ[3]*st.I3
    D = d.μ*st.I3
    return [S,E,I1,I2,I3,R,D]
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


    return SEI3RDynamics(α,β,γ,p,μ)
end

d=getParams(5.0,6.0,0.5,0.1,0.1,0.15,0.05,0.1,6.0,8.0,2.0)

state = [1000.0,0.0,0.0,0.0,0.0,0.0,0.0]

ntime = 300
nstates = Matrix{Float64}
for i in 1:300
  n = oneDay(n,d)
end
