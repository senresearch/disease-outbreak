using LsqFit
using CSV
include("../src/sirx.jl")
include("../src/evolution.jl")
# Hubei
hubeiFile = joinpath(@__DIR__, "..", "data",
     "covid19_china","time_series_covid19_confirmed_hubei.csv")
shandongFile = joinpath(@__DIR__, "..", "data",
          "covid19_china","time_series_covid19_confirmed_shandong.csv")

hubei = CSV.read(hubeiFile)
shandong = CSV.read(shandongFile)
hubeiPop = 57.0e6
shandongPop = 94.2e6

d = getParams(0.038,0.073,6.2,8.0,SIRX())
s0 = initialize(hubeiPop,d)
(s,ds) = evolve(hubeiPop.N,s0,d,100);
plotEvolution(hubeiPop.N,s,d,log10=true)

hubeiC = convert(Vector{Float64},hubei.ConfirmedCases[1:23])
shandongC = convert(Vector{Float64},shandong.ConfirmedCases[1:23])

hubeiFit = fitCaseModel(23,hubeiC,hubeiPop,6.2,8.0,[0.1,0.1,2.0])
shandongFit = fitCaseModel(23,shandongC,shandongPop,6.2,8.0,[0.5,0.045,15.0])
summary(hubeiFit)
summary(shandongFit)


plot(hubeiPop.N*s'[:,4],yaxis=:log,xaxis=:log,
     ylim=[10, 10^5],label="fitted")
plot!((hubei.ConfirmedCases),xaxis=:log,yaxis=:log,label="actual")
plot!(hubeiC,yaxis=:log,label="actual")
plot!(fitted(hubeiFit),yaxis=:log, label="estimated")
# savefig("hubei.pdf")

plot(shandongC,yaxis=:log,label="actual")
plot!(fitted(shandongFit),yaxis=:log, label="estimated")
