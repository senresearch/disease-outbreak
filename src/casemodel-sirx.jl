########################################################
# Functions related to fitting case series to SIRX model
########################################################

import Base.summary
using DataFrames

struct CaseModelFitResult
    κ::Float64
    κ0::Float64
    IXRatio::Float64
    inputs::NamedTuple
    fit::LsqFit.LsqFitResult
end

"""
fitCaseModel(nt::Int64,cases::Vector{Float64},
    N::Float64,R0Free::Float64,TInfected::Float64,
    p0::Vector{Float64})

- nt: Integer with length of series
- cases: Vector of case number series
- N: Population size
- R0Free: User supplied R0 in freely spreading population
- TInfected: Average time of infection
- p0: Floating point vector with initial values of estimates
"""
function fitCaseModel(nt::Int64,cases::Vector{Float64},
    N::Float64,R0Free::Float64,TInfected::Float64,
    p0::Vector{Float64})
    inputs = (nt=nt,C0=cases[1],N=N,R0Free=R0Free,TInfected=TInfected,C=cases)
    model(t,p) = logCaseModelUnknown(nt,exp.(p),N,cases[1],R0Free,TInfected)
    fit = curve_fit(model,(1:nt)*1.0,log.(cases),log.(p0))
    return CaseModelFitResult( exp(fit.param[1]), exp(fit.param[2]),
    exp(fit.param[3]), inputs, fit )
end

#######################################################
function fitCaseModel_ps(nt::Int64,cases::Vector{Float64},
    N::Float64,R0Free::Float64,TInfected::Float64,p0::Vector{Float64})
    inputs = (nt=nt,C0=cases[1],N=N,R0Free=R0Free,TInfected=TInfected,C=cases)
    logcases = log.(cases)
    f(p) = mse(logcases,logCaseModelUnknown(nt,exp.(p),N,cases[1],
               R0Free,TInfected),nt-3.0)
    fit = optimize(f,log.(p0),ParticleSwarm(lower=-3.0.*[1.0,1.0,0.0],
                   upper=3.0.*[0.0,0.0,1.0],n_particles=10),
                   Optim.Options(iterations=10000))
    # fit = optimize(f,log.(p0),NelderMead())
#    fit = optimize(f,log.(p0),-# 3.0.*[1.0,1.0,1.0],3.0.*[1.0,1.0,1.0],SAMIN())
    param = Optim.minimizer(fit)

    return ( κ=exp(param[1]), κ0= exp(param[2]),
             IXRatio=exp(param[3]), inputs=inputs, fit=fit )
end

function fitCaseModel_nm(nt::Int64,cases::Vector{Float64},
    N::Float64,R0Free::Float64,TInfected::Float64,p0::Vector{Float64})
    inputs = (nt=nt,C0=cases[1],N=N,R0Free=R0Free,TInfected=TInfected,C=cases)
    logcases = log.(cases)
    f(p) = mse(logcases,logCaseModelUnknown(nt,exp.(p),N,cases[1],
               R0Free,TInfected),nt-3.0)
    fit = optimize(f,log.(p0),NelderMead())

    param = Optim.minimizer(fit)

    return ( κ=exp(param[1]), κ0= exp(param[2]),
             IXRatio=exp(param[3]), inputs=inputs, fit=fit )
end

function mse(y::Vector{Float64},yhat::Vector{Float64},df::Float64)
    return sum( (y-yhat).^2  )/df
end


#######################################################

function summary(fit::CaseModelFitResult)
    r0eff = R0Eff(getParams(fit.κ,fit.κ0,fit.inputs.R0Free,
            fit.inputs.TInfected,SIRX()))
    df = exp.(DataFrame(LsqFit.confidence_interval(fit.fit)))
    rename!(df,[:lo,:hi])
    param = [ "κ", "κ0", "IXRatio" ]
    est = [ fit.κ, fit.κ0, fit.IXRatio ]
    insertcols!(df,1,parameter=param)
    insertcols!(df,3,est=est)
    return (coef=df,R0Eff=r0eff,mse=LsqFit.mse(fit.fit))
end

function fitted(fit::CaseModelFitResult,nt::Int64=0)
    if (nt==0)
        nt = fit.inputs.nt
    end
    return exp.(logCaseModelUnknown(nt,
                exp.(fit.fit.param), fit.inputs.N, fit.inputs.C0,
                fit.inputs.R0Free,fit.inputs.TInfected))
end

function estimatedStates(fit::CaseModelFitResult,nt::Int64=0)
    if (nt==0)
        nt = fit.inputs.nt
    end
    d = getParams(exp(fit.fit.param[1]),
                  exp(fit.fit.param[2]), fit.inputs.R0Free,
                  fit.inputs.TInfected, SIRX() )
    return estimatedStates( nt, fit.inputs.N, fit.inputs.C0,
           exp(fit.fit.param[3]), d )
end
