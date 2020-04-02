# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.1
#   kernelspec:
#     display_name: Julia 1.4.0
#     language: julia
#     name: julia-1.4
# ---

using Pkg

using Plots
using GraphRecipes, LightGraphs, Distributions, Distances, DataFrames, StatsBase

# +
function makeCitizens(n::Int64)
    x = rand(Uniform(),n);
    y = rand(Uniform(),n);
    citizenPos = [x y];
    citizenDists=pairwise(Euclidean(),citizenPos,dims=1)
    citizenStatus = DataFrame(ind=1:n,infected=zeros(Int64,n))
    return citizenStatus,citizenDists,citizenPos
end

function infectCitizen!(citizenStatus::DataFrame)
    citizenStatus[1,2]=1
end

# +
function updateStatus!(citizens::DataFrame,distances::Matrix{Float64},
        time::Int64,infectionDist::Float64=0.2,
        infectionWindow::Int64=1,infectionProb::Float64=0.2)
    whoCanInfect = citizens.ind[(citizens.infected.>0).&((time.-citizens.infected).<=infectionWindow)]
    whoSusceptible = citizens.ind[citizens.infected.==0]
    transmissionPossible = (distances[whoCanInfect,whoSusceptible].<=infectionDist)
    prob = 1.0.-(1.0.-infectionProb).^(sum(transmissionPossible,dims=1))
    # citizens.infected[whoSusceptible[rand.(Bernoulli.(prob))]] = time
    citizens.infected[citizens.ind[whoSusceptible[rand.(Bernoulli.(prob))[1,:]]]] .= time
end

function f(x)
    if(x==true)
        return 0 
    else
        return 1
    end
end
   
# -

function followCitizens!(citizens::DataFrame,distances::Matrix{Float64},ntime::Int64=300,
                        infectionDist::Float64=0.1,
                        infectionWindow::Int64=1,infectionProb::Float64=0.2)
    print(round(size(citizens,1)*pi*infectionDist^2*infectionWindow*infectionProb,digits=2))
    for i in 2:ntime
         updateStatus!(citizens,distances,i,infectionDist,
            infectionWindow,infectionProb)
    end
end

function uninfected(status::DataFrame) 
    return mean(status.infected.>0)
end


# +
function newInfections(status::DataFrame)
    histogram(status.infected,bins=range(1,maximum(status.infected),step=1))
end

function cumulativeInfections(status::DataFrame)
    plot(1:365,cumsum(counts(status.infected,1:365))./size(status,1),ylim=(0,1))
end

function infectionSpread(status::DataFrame,pos::Matrix{Float64})
    spread = @animate for i in 1:365
      scatter(pos[:,1],pos[:,2],color=f.((status.infected.<=i).&(status.infected.>0)))
    end
    return spread
end
# -

n=10000
status,d,pos = makeCitizens(n);
infectCitizen!(status)
followCitizens!(status,d,365,0.05,14,0.0025)

uninfected(status)

newInfections(status)

cumulativeInfections(status)

gif(infectionSpread(status,pos),fps=15)

n=10000
status1,d,pos = makeCitizens(n);
infectCitizen!(status1)
followCitizens!(status1,d,365,0.05,14,0.002)
uninfected(status1)

newInfections(status1)

cumulativeInfections(status1)

gif(infectionSpread(status1,pos),fps=15)
