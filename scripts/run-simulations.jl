using DiseaseOutbreak
using Plots

#####################################
# more infectious
#####################################
n=10000
status,d,pos = makeCitizens(n);
#TODO: confirm with Saunak: infectCitizens needs second argument: nInfect
infectCitizens!(status, 10)
followCitizens!(status,d,365,0.05,14,0.0025)

uninfected(status)

newInfections(status)

cumulativeInfections(status)

gif(infectionSpread(status,pos),fps=15)

#####################################
# less infectious
#####################################
n=10000
status1,d1,pos1 = makeCitizens(n);
infectCitizens!(status1, 10)
followCitizens!(status1,d1,365,0.05,14,0.002)
uninfected(status1)

newInfections(status1)

cumulativeInfections(status1)

gif(infectionSpread(status1,pos1),fps=15)
