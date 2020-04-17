using CSV, DataFrames,StatsBase
include(joinpath(@__DIR__, "..", "src", "distancing.jl"))

memphis = CSV.read(joinpath(@__DIR__, "..", "data","Memphis_residences.csv"),type=Float64,
                    missingstring="NA")

names(memphis)
memphispos = convert(Matrix{Float64},dropmissing(memphis[:,3:4]))

n = size(memphispos,1);
idx = sample(1:n,20000)
sampleMemphispos = memphispos[idx,:]
status,d,pos = makeCitizens(sampleMemphispos)
infectCitizens!(status,10)
followCitizens!(status,d,100,0.5,14,0.000015)

uninfected(status)

newInfections(status)

cumulativeInfections(status)

gif(infectionSpread(status,pos,80),joinpath(@__DIR__, "..", "plots", "memphis.gif"),fps=15)
