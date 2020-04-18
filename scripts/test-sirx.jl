include("../src/sirx.jl")
include("../src/evolution.jl")

# Hubei
using CSV
f = joinpath(@__DIR__, "..", "data", "covid19_china/time_series_covid19_confirmed_hubei.csv")

hubei = CSV.read(f)
N=57.1e6
d = getParams(0.038,0.073,6.2,8.0,SIRX())
s0 = initialize(N,hubei.ConfirmedCases[1]*1.0,2.55,d)
(s,ds) = evolve(N,s0,d,100);
plotEvolution(N,s,d,log10=true)

plot(N*s'[:,4],yaxis=:log,xaxis=:log,ylim=[10, 10^5],label="fitted")
plot!((hubei.ConfirmedCases),xaxis=:log,yaxis=:log,label="actual")
# savefig("hubei.pdf")
