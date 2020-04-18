include("../src/sirx.jl")
include("../src/evolution.jl")

# Hubei
using CSV
f = joinpath(@__DIR__, "..", "data",
     "covid19_china/time_series_covid19_confirmed_hubei.csv")

hubei = CSV.read(f)
N=57.1e6
d = getParams(0.038,0.073,6.2,8.0,SIRX())
s0 = initialize(N,hubei.ConfirmedCases[1]*1.0,2.55,d)
(s,ds) = evolve(N,s0,d,100);
plotEvolution(N,s,d,log10=true)

plot(N*s'[:,4],yaxis=:log,xaxis=:log,ylim=[10, 10^5],label="fitted")
plot!((hubei.ConfirmedCases),xaxis=:log,yaxis=:log,label="actual")
# savefig("hubei.pdf")
function caseModel(t::Vector{Float64},params::Vector{Float64},
    N::Float64,C0::Float64,IXRatio::Float64,R0::Float64,TI::Float64)
    d = getParams(params[1],params[2],R0,TI,SIRX())
    s0 = initialize(N,C0,IXRatio,d)
    (d,ds) =evolve(N,s0,d,length(t))
    return N.*d[4,:]
end

z = caseModel(convert(Vector{Float64},
              (1:length(hubei.ConfirmedCases))),
              [0.038,0.073],N,
              hubei.ConfirmedCases[1]*1.0,
              2.55,6.2,8.0)

plot(N.*z,yaxis=:log,xaxis=:log,ylim=[10, 10^5],label="fitted")
plot!((hubei.ConfirmedCases),xaxis=:log,yaxis=:log,label="actual")

function model(x::Vector{Float64},p::Vector{Float64})
    return log.(caseModel(x,exp.(p),N,
                hubei.ConfirmedCases[1]*1.0,2.55,6.2,8.0))
end

using LsqFit

x = convert(Vector{Float64},(1:length(hubei.ConfirmedCases)))
y = log.(hubei.ConfirmedCases)
fit = curve_fit(model,x,y,[log(0.038),log(0.073)])

fit
exp.(fit.param)
