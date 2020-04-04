using CSV, DataFrames,StatsBase
include("src/distancing.jl")

memphis = CSV.read("data/Memphis_residences.csv",
                   type=Float64,missingstring="NA")

names(memphis)
memphispos = convert(Matrix{Float64},dropmissing(memphis[:,3:4]))

n = size(memphispos,1);
idx = sample(1:n,10000)
sampleMemphispos = memphispos[idx,:]
status,d,pos = makeCitizens(sampleMemphispos)
infectCitizens!(status,5)
followCitizens!(status,d,365,0.1,14,0.00005)

uninfected(status)

newInfections(status)

cumulativeInfections(status)

gif(infectionSpread(status,pos,360),fps=15)
