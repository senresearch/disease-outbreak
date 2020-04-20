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

hubeiC = convert(Vector{Float64},hubei.ConfirmedCases[1:23])
shandongC = convert(Vector{Float64},shandong.ConfirmedCases[1:23])

hubeiFit = fitCaseModel(23,hubeiC,hubeiPop,
                        6.2,8.0,[0.1,0.1,2.0])
shandongFit = fitCaseModel(23,shandongC,shandongPop,
                           6.2,8.0,[0.5,0.045,15.0])
summary(hubeiFit)
summary(shandongFit)

plot(hubeiC,xaxis=:log,yaxis=:log,label="actual")
plot!(caseModel(23,hubeiPop,hubeiC[1],2.26,
                getParams(0.000,0.084,6.2,8.0,SIRX())),
                yaxis=:log,label="paper")
plot!(fitted(hubeiFit),yaxis=:log, label="estimated")
# savefig("hubei.pdf")

plot(shandongC,yaxis=:log,label="actual")
plot!(caseModel(23,shandongPop,shandongC[1],9.66,
                getParams(0.309,0.042,6.2,8.0,SIRX())),
                yaxis=:log,label="paper")
plot!(fitted(shandongFit),yaxis=:log, label="estimated")
