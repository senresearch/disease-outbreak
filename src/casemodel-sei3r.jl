using Optim
using Calculus
import DiseaseOutbreak.fitCaseModel

logit(x::Float64) = log(x/(1-x))
invlogit(x::Float64) = exp(x)/(1+exp(x))

function poissonLogLik( x::Vector{Int64}, λ::Vector{Float64})
    return sum(poissonLogLik.(x,λ))
end

function poissonLogLik(x::Int64,λ::Float64)
    return logpdf(Poisson(λ),x)
end

"""
function fitCaseModel(cases::Vector{Int64},N::Int64,
    ρ::Float64,ϕ0::Float64,E::Float64,d0::SEI3R)
"""
function fitCaseModel(cases::Vector{Int64},N::Int64,
    ρ::Float64,ϕ0::Float64,E::Float64,d0::SEI3R)

    d = deepcopy(d0)

    ntime = length(cases)
    # dvec = fill(d0,ntime)

    function f(param::Vector{Float64})
        d.β[1] = invlogit(param[2] + logit(d0.β[1]) )
        s = estimatedStates(ntime,N,cases[1],exp(param[3]),
                            invlogit(param[1]),d)
        caseIntensity = s.I1*invlogit(param[1])
        # println(caseIntensity)
        return -poissonLogLik(diff(cases),caseIntensity[2:end])
        # return caseIntensity
    end

    fit = optimize(f,[logit(ρ),logit(ϕ0),log(E)],NelderMead())
    d.β[1] = invlogit(fit.minimizer[2] + logit(d0.β[1]) )
    info = Calculus.hessian(f,fit.minimizer)
    return (cases=cases,N=N,ρ=invlogit(fit.minimizer[1]),
            E=exp(fit.minimizer[3]),d=d,info=info,fit=fit)
end

###########
"""
function fitCaseModel(cases::Vector{Int64},x::Matrix{Float64},
    N::Int64,ρ::Float64,ϕ0::Float64,E::Float64,d0::SEI3R)

Fit a model to confirmed cases where the β1 transmission parameter is allowed to depend on covariates x.
"""
function fitCaseModel(cases::Vector{Int64},x::Matrix{Float64},
    N::Int64,ρ::Float64,ϕ0::Float64,E::Float64,d0::Vector{SEI3R})

    d = deepcopy(d0)

    ntime = length(cases)
    # dvec = fill(d0,ntime)

    function f(param::Vector{Float64})
        d.β[1] = invlogit(param[2] + logit(d0.β[1]) )
        s = estimatedStates(ntime,N,cases[1],exp(param[3]),
                            invlogit(param[1]),d)
        caseIntensity = s.I1*invlogit(param[1])
        # println(caseIntensity)
        return -poissonLogLik(diff(cases),caseIntensity[2:end])
        # return caseIntensity
    end

    fit = optimize(f,[logit(ϕ),logit(ϕ1),log(E)],NelderMead())
    d.β[1] = invlogit(fit.minimizer[2] + logit(d0.β[1]))
    info = Calculus.hessian(f,fit.minimizer)
    return (cases=cases,N=N,ϕ=invlogit(fit.minimizer[1]),
            E=exp(fit.minimizer[3]),d=d,info=info,fit=fit)
end


function predictCases(fit,ntime)
    s = estimatedStates(ntime,fit.N,fit.cases[1],fit.E,fit.ϕ,fit.d)
    return cumsum(s.I1*fit.ϕ)
end
