include("../src/sirx.jl")
include("../src/evolution.jl")

# Hubei
using CSV
hubeiFile = joinpath(@__DIR__, "..", "data",
     "covid19_china","time_series_covid19_confirmed_hubei.csv")

hubei = CSV.read(hubeiFile)
hubeiPop=57.1e6
d = getParams(0.038,0.073,6.2,8.0,SIRX())
s0 = initialize(hubeiPop,hubei.ConfirmedCases[1]*1.0,2.55,d)
(s,ds) = evolve(hubeiPop,s0,d,100);
plotEvolution(hubeiPop,s,d,log10=true)

using LsqFit
using Optim

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

plot(hubeiPop*s'[:,4],yaxis=:log,xaxis=:log,
     ylim=[10, 10^5],label="fitted")
plot!((hubei.ConfirmedCases),xaxis=:log,yaxis=:log,label="actual")
plot!(caseModel(x,exp.(fit.param),hubeiPop,
                       hubei.ConfirmedCases[1]*1.0,
                       6.2,8.0),yaxis=:log, label="estimated")
# savefig("hubei.pdf")

g(p) = ss(y,x,p)
git = optimize(g,log.([0.1,0.1,1.0]),NelderMead())

exp.(git.minimizer)
