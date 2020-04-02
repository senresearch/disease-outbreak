include("distancing.jl")

#####################################
# more infectious
#####################################
n=10000
status,d,pos = makeCitizens(n);
infectCitizen!(status)
followCitizens!(status,d,365,0.05,14,0.0025)

uninfected(status)

newInfections(status)

cumulativeInfections(status)

gif(infectionSpread(status,pos),fps=15)

#####################################
# less infectious
#####################################
n=10000
status1,d,pos = makeCitizens(n);
infectCitizen!(status1)
followCitizens!(status1,d,365,0.05,14,0.002)
uninfected(status1)

newInfections(status1)

cumulativeInfections(status1)

gif(infectionSpread(status1,pos),fps=15)
