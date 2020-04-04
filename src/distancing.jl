## load packages
using Plots
using Distributions, Distances, DataFrames, StatsBase

```
makeCitizens(n::Int64)

Create a population of `n` individuals uniformly distributed in a rectangle.

Returns the spatial positions of the citizens, their Euclidean distances from
each other, and initializes their disease status to being unaffected (0).

This function can be generalized to include other spatial distributions, and
other ways to calculate distance between citizens.
```
function makeCitizens(n::Int64)
    x = rand(Uniform(),n);
    y = rand(Uniform(),n);
    citizenPos = [x y];
    citizenDists=pairwise(Euclidean(),citizenPos,dims=1)
    citizenStatus = DataFrame(ind=1:n,infected=zeros(Int64,n))
    return citizenStatus,citizenDists,citizenPos
end

```makeCitizens(citizenPos::Matrix{Float64})
```
function makeCitizens(citizenPos::Matrix{Float64})
    n = size(citizenPos,1)
    citizenDists = pairwise(Euclidean(),citizenPos,dims=1)
    citizenStatus = DataFrame(ind=1:n,infected=zeros(Int64,n))
    return citizenStatus,citizenDists,citizenPos
end
```
infectCitizens!(citizenStatus::DataFrame)

Picks `nInfect` random individuals in `citizenStatus` to be infected.  This
function could be generalized by picking someone in a defined region or
with defined characteristics.
```
function infectCitizens!(citizenStatus::DataFrame,nInfect::Int64)
    n = size(citizenStatus,1)
    citizenStatus[sample(1:n,nInfect),2] .= 1
end

```
updateStatus!(status::DataFrame,distances::Matrix{Float64},
             time::Int64,infectionDist::Float64,infectionWindow::Int64,
             infectionProb::Float64)

Updates the infection `status` data frame at a given `time` using information
on `infectionDist` (radius of neighborhood which a person can infect),
`infectionWindow` (how long an infected person is infectious), `infectionProb`
(the probability of infecting a susceptible person in a given time period).

The infection model can be tweaked to make it more realistic.
```
function updateStatus!(status::DataFrame,distances::Matrix{Float64},
        time::Int64,infectionDist::Float64,
        infectionWindow::Int64,infectionProb::Float64)

    # individuals who are still infectious
    whoCanInfect = status.ind[(status.infected.>0).&
                   ((time.-status.infected).<=infectionWindow)]
    # individuals who are susceptible
    whoSusceptible = status.ind[status.infected.==0]
    # which individuals can infect whom
    transmissionPossible = (distances[whoCanInfect,whoSusceptible].<=infectionDist)
    # calculate probability of infection for each person in a given time period
    prob = 1.0.-(1.0.-infectionProb).^(sum(transmissionPossible,dims=1))
    # randomly decide based on infection probability, who is newly infected
    status.infected[status.ind[whoSusceptible[rand.(Bernoulli.(prob))[1,:]]]] .= time
end

```
followCitizens!(citizens::DataFrame,distances::Matrix{Float64},ntime::Int64=365,
                        infectionDist::Float64,
                        infectionWindow::Int64,infectionProb::Float64)

Follow the citizens in the data drame `citizens` for `ntime` time periods assuming
`infectionDist`, `infectionWindow` an `infectionProb` as in updateStatus function.
```
function followCitizens!(citizens::DataFrame,distances::Matrix{Float64},ntime::Int64=365,
                        infectionDist::Float64=0.1,
                        infectionWindow::Int64=1,infectionProb::Float64=0.2)
    # print expected number of infections per infected individual
    print(round(size(citizens,1)*pi*infectionDist^2*infectionWindow*infectionProb,digits=2))
    for i in 2:ntime
         updateStatus!(citizens,distances,i,infectionDist,
            infectionWindow,infectionProb)
    end
end

###
## convert trues to 1
###
function trueone(x)
    if x==true
        return 0
    else
        return 1
    end
end

###################################################


```
uninfected(status::DataFrame)

Returns number uninfected at end of follow up period.
```
function uninfected(status::DataFrame)
    return mean(status.infected.>0)
end

```
newInfections(status::DataFrame)

Make histogram of new infections.
```
function newInfections(status::DataFrame)
    histogram(status.infected,bins=range(1,maximum(status.infected),step=1))
end

```
cumulativeInfections(status::DataFrame)

Plot cumulative number of infection as a proportion of the population.
```
function cumulativeInfections(status::DataFrame)
    plot(1:365,cumsum(counts(status.infected,1:365))./size(status,1),ylim=(0,1))
end

```
infectionSpread(status::DataFrame,pos::Matrix{Float64})

Make a animation of the spread of the infection using the `status` data and
spatial positions in `pos`.
```
function infectionSpread(status::DataFrame,pos::Matrix{Float64},ntime::Int64=365)
    spread = @animate for i in 1:ntime
      scatter(pos[:,1],pos[:,2],
              color=trueone.((status.infected.<=i).&(status.infected.>0)))
    end
    return spread
end
