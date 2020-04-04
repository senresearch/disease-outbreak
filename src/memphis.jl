using CSV, DataFrames
include("distancing.jl")

memphis = CSV.read("data/Memphis_residences.csv",type=Float64)

names(memphis)
memphispos = convert(Matrix{Float64},dropmissing(memphis[:,3:4]))

n = size(memphispos,1);

status,d,pos = makeCitizens(memphispos)
infectCitizen!(status)
followCitizens!(status,d,365,0.05,14,0.0025)

uninfected(status)

newInfections(status)

cumulativeInfections(status)

gif(infectionSpread(status,pos),fps=15)
