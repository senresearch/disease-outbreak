########################################################
# Functions related to fitting case series to SIRX model
########################################################

include("sirx.jl")

struct CaseModelFitResult
    κ::Float64
    κ0::Float64
    IXRatio::Float64
    inputs::NamedTuple
    fit::LsqFit.LsqFitResult
end

function fitCaseModel(nt::Int64,cases::Vector{Float64},
    N::Float64,R0Free::Float64,TInfected::Float64,
    p0::Vector{Float64})
    inputs = (nt=nt,C0=cases[1],N=N,R0Free=R0Free,TInfected=TInfected,C=cases)
    model(t,p) = logCaseModelUnknown(nt,exp.(p),N,cases[1],R0Free,TInfected)
    fit = curve_fit(model,(1:nt)*1.0,log.(cases),log.(p0))
    return CaseModelFitResult( exp(fit.param[1]), exp(fit.param[2]),
    exp(fit.param[3]), inputs, fit )
end


function summary(fit::CaseModelFitResult)
    r0eff = R0Eff(getParams(fit.κ,fit.κ0,fit.inputs.R0Free,
            fit.inputs.TInfected,SIRX()))
    return (κ=fit.κ,κ0=fit.κ0,IXRatio=fit.IXRatio,R0Eff=r0eff)
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


function plotfit(fit::CaseModelFitResult,nt::Int64=0)
    if (nt==0)
        nt = fit.inputs.nt
    end
    plot(fit.inputs.C,yaxis=:log,seriestype=:scatter,
                color=:black,label="actual")
    plot!(estimatedStates(fit,nt)[:I],
                yaxis=:log,label="infected")
    plot!(fitted(fit,nt),yaxis=:log, label="fitted",color=:blue)
end
