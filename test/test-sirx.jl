using CSV, Plots
using DiseaseOutbreak

# Hubei and Shandong
hubeiFile = joinpath(@__DIR__, "..", "data",
     "covid19_china","time_series_covid19_confirmed_hubei.csv")
shandongFile = joinpath(@__DIR__, "..", "data",
          "covid19_china","time_series_covid19_confirmed_shandong.csv")

hubei = CSV.read(hubeiFile)
shandong = CSV.read(shandongFile)
hubeiPop = 57.0e6
shandongPop = 94.2e6

hubeiC = convert(Vector{Float64},hubei.ConfirmedCases[1:22])
hubeiFit = fitCaseModel(22,hubeiC,hubeiPop,
                        6.2,8.0,[0.1,0.1,2.0])
summary(hubeiFit)
plotFit(hubeiFit,40)
# savefig("hubei.pdf")
pyplotFit(hubeiFit, 40, log10="semilog", grid=true, fsize=[12, 8])
pyplotFit(hubeiFit, 40, log10="loglog", grid=true, fsize=[12, 8])
pyplotFit(hubeiFit, 40, log10="null", grid=true, fsize=[12, 8])

shandongC = convert(Vector{Float64},shandong.ConfirmedCases[1:30])
shandongFit = fitCaseModel(30,shandongC,shandongPop,
                           6.2,8.0,[0.5,0.045,15.0])
summary(shandongFit)

plotFit(shandongFit,40)
plot!(estimatedStates(30,shandongPop,shandongC[1],9.66,
                getParams(0.309,0.042,6.2,8.0,SIRX()))[!,:X],
                yaxis=:log,label="paper")

pyplotFit(shandongFit, 40, log10="semilog", grid=true, fsize=[12, 8])
paperRef = estimatedStates(30,shandongPop,shandongC[1],9.66,
                getParams(0.309,0.042,6.2,8.0,SIRX()))[!,:X]

pyplotFit(shandongFit, 30, ref=paperRef, log10="semilog", grid=true,
          fsize=[12, 8])

fitted(shandongFit,40)
