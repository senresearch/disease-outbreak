using LsqFit
using CSV
include("../src/sirx.jl")
include("../src/evolution.jl")
# Hubei
hubeiFile = joinpath(@__DIR__, "..", "data",
     "covid19_china","time_series_covid19_confirmed_hubei.csv")

hubei = CSV.read(hubeiFile)
d = getParams(0.038,0.073,6.2,8.0,SIRX())
hubeiPop = SIRXPopulation(57.0e6,hubei.ConfirmedCases[1]*1.0,2.55)
s0 = initialize(hubeiPop,d)
(s,ds) = evolve(hubeiPop.N,s0,d,100);
plotEvolution(hubeiPop.N,s,d,log10=true)

t = convert(Vector{Float64},(1:length(hubei.ConfirmedCases)))
y = convert(Vector{Float64},hubei.ConfirmedCases)

hubeiPop.N
fit = fitCaseModel(t,y,hubeiPop.N,6.2,8.0,[0.1,0.1,2.0])
summary(fit)

plot(hubeiPop.N*s'[:,4],yaxis=:log,xaxis=:log,
     ylim=[10, 10^5],label="fitted")
plot!((hubei.ConfirmedCases),xaxis=:log,yaxis=:log,label="actual")
plot!(caseModel(x,exp.(fit.param),hubeiPop,
                       hubei.ConfirmedCases[1]*1.0,
                       6.2,8.0),yaxis=:log, label="estimated")
# savefig("hubei.pdf")
