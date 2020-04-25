using Optim
import DiseaseOutbreak.fitCaseModel

logit(x::Float64) = log(x/(1-x))
invlogit(x::Float64) = exp(x)/(1+exp(x))

function poissonLogLik( x::Vector{Int64}, λ::Vector{Float64})
    return sum(poissonLogLik.(x,λ))
end

function poissonLogLik(x::Int64,λ::Float64)
    return logpdf(Poisson(λ),x)
end

function fitCaseModel(cases::Vector{Int64},N::Int64,
    ϕ::Float64,ϕ1::Float64,d0::SEI3R)

    d = deepcopy(d0)

    ntime = length(cases)
    # dvec = fill(d0,ntime)

    function f(param::Vector{Float64})
        d.β[1] = invlogit(param[2] + logit(d0.β[1]) )
        s = estimatedStates(ntime,N,cases[1],invlogit(param[1]),d)
        caseIntensity = s.I1*invlogit(param[1])
        # println(caseIntensity)
        return -poissonLogLik(diff(cases),caseIntensity[2:end])
        # return caseIntensity
    end

    fit = optimize(f,logit.([ϕ,ϕ1]),NelderMead())
    d.β[1] = invlogit(fit.minimizer[2] + logit(d0.β[1]) )
    return (cases=cases,N=N,ϕ=invlogit(fit.minimizer[1]),d=d,fit=fit)
end

function predictCases(fit,ntime)
    s = estimatedStates(ntime,fit.N,fit.cases[1],fit.ϕ,fit.d)
    return cumsum(s.I1*fit.ϕ)
end
