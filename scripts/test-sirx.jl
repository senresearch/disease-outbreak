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
