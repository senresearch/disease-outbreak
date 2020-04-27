# include("sei3r.jl")
# include("evolution.jl")

using XLSX
using DiseaseOutbreak
# TODO: make numerical tests to replace plot tests
using Test

# default parameters for RCAN model
d0 = getParams(5.0,6.0,0.5,0.1,0.1,0.15,0.06,6.0,8.0,0.4,SEI3R())
# social distancing reducing transmission probs by 1/2
d1 = getParams(5.0,6.0,0.5/2,0.1/2,0.1/2,0.15,0.06,6.0,8.0,0.4,SEI3R())

N = 1000.0
state0 = initialize(1.0/N,d0)
(states0,deltaStates0) = evolve(N,state0,d0,300)
plotEvolution(N,states0,d0)
# savefig("dist300.pdf")
state0 = initialize(1.0/N,d0)
(states0,deltaStates0) = evolve(N,state0,d0,60)
(states1,deltaStates1) = evolve(N,states0[:,60],d1,241)
plotEvolution(N,[states0 states1[:,2:end]],d0)
# savefig("dist60.pdf")
state0 = initialize(1.0/N,d0)
(states0,deltaStates0) = evolve(N,state0,d0,30)
(states1,deltaStates1) = evolve(N,states0[:,30],d1,271)
plotEvolution(N,[states0 states1[:,2:end]],d0)
# savefig("dist30.pdf")
state1 = initialize(1.0/N,d1)
(states1,deltaStates1) = evolve(N,state1,d1,300)
plotEvolution(N,states1,d1)
# savefig("dist0.pdf")


sheet = XLSX.readdata(joinpath("../memphis-covid-data",
             "covid_status_updates/",
             "COVID Status April 23.xlsx"),
             "Status","A1:C49")
using CSV
enc = CSV.read(joinpath("../memphis-covid-data",
               "mobility_data","unacast-shelby.csv"))
memCases = convert(Vector{Int64},sheet[2:49,2])
memPop = 1400000

fit = fitCaseModel(memCases[13:end],memPop,0.2,0.5,10.0,d0)
plot(memCases[13:end],legend=:bottomright,label="cumulative cases")
prd = predictCases(fit,60)
plot!(prd,label="fit")

length(memCases[13:end])
encMetric = reverse(enc.encountersMetric[1:36])
x = [ones(36) encMetric]
fit1 = fitCaseModel((memCases[13:end]),x[1:35,:],memPop,0.3,10.0,d0)
plot(memCases[13:end],legend=:bottomright,label="cumulative cases")
prd = predictCases(fit,60)
plot!(prd,label="fit0")

prd1 = predictCases(fit1,x[1:35,:])
plot!(prd1,label="fit1")

x2 = ones(35,2);
x2[1:13,2] = zeros(13)
fit2 = fitCaseModel(memCases[13:end],x2,memPop,0.01,10.0,d0)
prd2 = predictCases(fit2,x2)
plot!(prd2,label="fit2")

############

tn = CSV.read("/Users/sen/covid-memphis/tn.csv")

tnCases = reverse(tn.ConfirmedCases[1:42])
tnPop = 6800000
fittn = fitCaseModel(tnCases,tnPop,0.01,0.5,100.0,d0)
plot(tnCases,legend=:bottomright,yaxis=:log,label="cumulative cases")
prdtn = predictCases(fittn,60)
plot!(prdtn,label="fit0",yaxis=:log)
xtn = ones(length(tnCases)-1,2)
xtn[1:10,2] = zeros(10)

fittn1 = fitCaseModel(tnCases,xtn,tnPop,0.02,1000.0,d0)
xtn1 = ones(200,2)
xtn1[1:10,2] = zeros(10)
prdtn1 = predictCases(fittn1,xtn1)
plot!(prdtn1,label="fit2",yaxis=:log)

plot(diff(tnCases),legend=:bottomright)
