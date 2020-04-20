# test-sirx.jl is a script that test the model SIRX
include(realpath(joinpath((@__DIR__), "..", "src", "sirx.jl")))
include(realpath(joinpath((@__DIR__), "..", "src", "evolution.jl")))
include(realpath(joinpath((@__DIR__), "..", "src", "sirPlot.jl")))

using LsqFit
using CSV, Plots


# Read dataset
hubeiFile = joinpath(@__DIR__, "..", "data", "covid19_china",
                    "time_series_covid19_confirmed_hubei.csv")
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


Plots.plot(hubeiC,xaxis=:log,yaxis=:log,label="actual")
Plots.plot!(caseModel(23,hubeiPop,hubeiC[1],2.26,
                getParams(0.000,0.084,6.2,8.0,SIRX())),
                yaxis=:log,label="paper")
Plots.plot!(fitted(hubeiFit),yaxis=:log, label="estimated")

# savefig("hubei.pdf")

plot(shandongC,yaxis=:log,label="actual")
plot!(caseModel(23,shandongPop,shandongC[1],9.66,
                getParams(0.309,0.042,6.2,8.0,SIRX())),
                yaxis=:log,label="paper")
plot!(fitted(shandongFit),yaxis=:log, label="estimated")



# BASED ON THE LAST SCRIPT

# Set initial parameters
# getParams(κ, κ0, α, R0Free, TInfected, SIRX)
d = getParams(0.038, 0.073, 6.2, 8.0, SIRX())
hubeiPop = SIRXPopulation(57.0e6, hubei.ConfirmedCases[1]*1.0, 2.55)
# Get initial states S-I-R-X
s0 = initialize(hubeiPop,hubeiC[1],2.26, d)
# Get new states after n = 100 days
s, ds = evolve(hubeiPop, s0, d, 100);

# This plot use the function in sirPlot.jl
# it is possible to save with the keyword (e.g., savefigure = "my_plot.pdf")
# it also possible to change the size of the figure with fsize = [10, 7]
myPlotEvolution(hubeiPop, s, d, log10 = "semilog", grid = true)


t = convert(Vector{Float64},(1:length(hubei.ConfirmedCases)))
y = convert(Vector{Float64},hubei.ConfirmedCases)


fit = fitCaseModel(t, y, hubeiPop.N, 6.2, 8.0, [0.1,0.1,2.0])
summary(fit)

hubeiFit = fitCaseModel(23, hubeiC, hubeiPop,
                        6.2, 8.0, [0.1,0.1,2.0])
summary(hubeiFit)

# This plot use the function in sirPlot.jl
# it is possible to save with the keyword (e.g., savefigure = "my_plot.pdf")
# it also possible to change the size of the figure with fsize = [10, 7]
plotEvoFitted(hubeiPop, s, d, hubeiC,
              log10 = "loglog", grid = true)


# savefig("hubei.pdf") saving function has been added to the plot function
