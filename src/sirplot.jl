##########################################################
# Functions related to plotting case series from SIR model
##########################################################

using PyCall, PyPlot, Plots
pygui(:tk)

include(realpath(joinpath((@__DIR__), "sirx.jl")))
include(realpath(joinpath((@__DIR__), "casemodelfit.jl")))


"""
plotEvolution(N::Float64,states::Matrix{Float64},
    d::Dynamics;log10::Bool=false)

Generate a plot of states evolved for the SIR model.

- N: the population size
- states: a matrix of states evoloved over time. See the function, 'evolve'.
- d: a type of Dynamics of the SIR model
- log10: a keyword argument (default is false) of scaling state values by log10

"""

function plotEvolution(N::Float64,states::Matrix{Float64},
    d::Dynamics;log10::Bool=false)
    if (log10==false)
        Plots.plot(N*Matrix(states'),label=stateNames(d))
    else
        Plots.plot(N*Matrix(states'),label=stateNames(d),ylab=:log10)
    end
end


"""
plotfit(fit::CaseModelFitResult,nt::Int64=0)

- fit: output from fitCaseModel
- nt: how long we want to extend predictions (use zero for not extrapolating)
"""
function plotFit(fit::CaseModelFitResult,nt::Int64=0)
    if (nt==0)
        nt = fit.inputs.nt
    end
    Plots.plot(fit.inputs.C,yaxis=:log,seriestype=:scatter,
                color=:black,label="actual")
    Plots.plot!(estimatedStates(fit,nt)[!,:I],
                yaxis=:log,label="infected")
    Plots.plot!(fitted(fit,nt),yaxis=:log, label="fitted",color=:blue)
end


"""
pyplotEvolution(N::Float64, states::Matrix{Float64},
              d::Dynamics;  log10 = "Null", figTitle = "",
              grid = "off")

Generate a plot of states evolved for the SIR model.

# Arguments
- N: the population size
- states: a matrix of states evoloved over time. See the function, `evolve`.
- d: a type of Dynamics of the SIR model
- log10: a keyword argument (default is "Null") where "semilog" makes a plot
          with log scaling on the y axis, "loglog" makes a plot with log scaling
          on the x axis and y axis, and "Null"makes a plot without log scaling.
- figTitle: a keyword argument (default is empty) that add a title to the figure.
- grid: a keyword argument (default is `false`) that adds a grid on the plotting if
         it is `true`.
- fsize: a keyword argument (default is [5, 5]]) that set the size of the figure.
- savefigure: a keyword argument (default is empty "") that save the figure in
               a given path.(e.g., savefigure = "../images/Plot.png") The extension
               of the filename indicates the format that is used to save the figure
               such as .pdf, .png, .jpg,...
"""

function pyPlotEvolution(N::Float64, states::Matrix{Float64},
                         d::Dynamics; log10 = "Null", figTitle = "",
                         grid = false, fsize= [10, 7], savefigure = "")

    fig, (ax) = PyPlot.subplots(nrows=1, figsize=(fsize[1], fsize[2]))

    ###############
    #  PLOT DATA  #
    ###############

    if log10 == "semilog"
        ax.semilogy(N*Matrix(states'))
    elseif log10 == "loglog"
        ax.loglog(N*Matrix(states'))
    elseif log10 == "Null"
        ax.plot(N*Matrix(states'))
    else
        println("Incorrect log10 option, no log is applied.")
        ax.plot(N*Matrix(states'))
    end


    #####################
    #  AXES AND LABELS  #
    #####################

    ax.spines["top"].set_visible(false) # Hide the top edge of the axis
    ax.spines["right"].set_visible(false) # Hide the right edge of the axis
    ax.xaxis.set_ticks_position("bottom")
    # ax1.spines["left"].set_position(("axes",-0.03)) # Offset the left scale from the axis
    # ax1.spines["bottom"].set_position(("axes",-0.05)) # Offset the bottom scale from the axis



    ax.set_xlabel("Time")
    ax.set_ylabel("Confirmed")

    # axis("tight")

    PyPlot.grid(grid)
    PyPlot.legend(stateNames(d)[:], loc="center right", fancybox="true")
    PyPlot.title(figTitle)

    if !isempty(savefigure)
        PyPlot.savefig(savefigure)
    end
    PyPlot.display_figs()

end



"""
pyplotFit(fit::CaseModelFitResult, nt::Int64 = 0;
        log10 = "Null", figTitle = "", grid = false,
        fsize= [5, 5], savefigure = "", )

Generate a plot of fitted state evolution for the SIR model.

# Arguments
- pop: a type SIRXPopulation that contains population size
- states: a matrix of states evoloved over time. See the function, `evolve`.
- d: a type of Dynamics of the SIR model
- actualState: actual case numbers
- log10: a keyword argument (default is "Null") where "semilog" makes a plot
          with log scaling on the y axis, "loglog" makes a plot with log scaling
          on the x axis and y axis, and "Null"makes a plot without log scaling.
- figTitle: a keyword argument (default is empty) that add a title to the figure.
- grid: a keyword argument (default is `false`) that adds a grid on the plotting if
         it is `true`.
- fsize: a keyword argument (default is [5, 5]]) that set the size of the figure.
- savefigure: a keyword argument (default is empty "") that save the figure in
               a given path.(e.g., savefigure = "../images/Plot.png") The extension
               of the filename indicates the format that is used to save the figure
               such as .pdf, .png, .jpg,...
"""

function pyplotFit(fit::CaseModelFitResult, nt::Int64 = 0;
                 log10 = "Null", figTitle = "", grid = false,
                 fsize= [10, 7], savefigure = "", ref = [0])


    fig, (ax) = PyPlot.subplots(nrows=1, figsize=(fsize[1], fsize[2]))

    ###############
    #  PLOT DATA  #
    ###############

    if (nt==0)
        nt = fit.inputs.nt
    end

    if log10 == "semilog"
        ax.semilogy(fit.inputs.C, linestyle ="-", color = "#A7A7A7",
                  marker="s", markersize = 15, markerfacecolor= "none",
                  label= "actual")
        ax.semilogy(estimatedStates(fit, nt)[!, :I], linestyle ="--", color = "red",
                  label= "infected")
        ax.semilogy(fitted(fit, nt), label="fitted")
        if (length(ref) > 1)
            ax.semilogy(ref, color = "green",  label="paper")
        end
        xlim(0, nt)

    elseif log10 == "loglog"
        ax.loglog(fit.inputs.C, linestyle ="-", color = "#A7A7A7",
                  marker="s", markersize = 15, markerfacecolor= "none",
                  label= "actual")
        ax.loglog(estimatedStates(fit, nt)[!, :I], linestyle ="--", color = "red",
                  label= "infected")
        ax.loglog(fitted(fit, nt), label="fitted")
        if (length(ref) > 1)
            ax.semilogy(ref, color = "green",  label="paper")
        end
        # xlim(0, log10(nt+2))

    elseif log10 == "null"
        ax.plot(fit.inputs.C, linestyle ="-", color = "#A7A7A7",
                  marker="s", markersize = 15, markerfacecolor= "none",
                  label= "actual")
        ax.plot(estimatedStates(fit, nt)[!, :I], linestyle ="--", color = "red",
                  label= "infected")
        ax.plot(fitted(fit, nt), label="fitted")
        if (length(ref) > 1)
            ax.semilogy(ref, color = "green",  label="paper")
        end
        xlim(0, nt)
    else
        println("Incorrect log10 option, no log is applied.")
        ax.plot(fit.inputs.C, linestyle ="-", color = "#A7A7A7",
                  marker="s", markersize = 15, markerfacecolor= "none",
                  label= "actual")
        ax.plot(estimatedStates(fit, nt)[!, :I], linestyle ="--", color = "red",
                  label= "infected")
        ax.plot(fitted(fit, nt), label="fitted")
        if (length(ref) > 1)
            ax.semilogy(ref, color = "green",  label="paper")
        end
        xlim(0, nt)
    end


    #####################
    #  AXES AND LABELS  #
    #####################

    ax.spines["top"].set_visible(false) # Hide the top edge of the axis
    ax.spines["right"].set_visible(false) # Hide the right edge of the axis
    ax.xaxis.set_ticks_position("bottom")
    # ax.spines["left"].set_position(("axes",-0.03)) # Offset the left scale from the axis
    # ax.spines["bottom"].set_position(("axes",-0.05)) # Offset the bottom scale from the axis
    ax.autoscale(false)


    ax.set_xlabel("Time")
    ax.set_ylabel("Confirmed")

    # axis("tight")

    PyPlot.grid(grid)
    PyPlot.legend( loc="lower center")
    PyPlot.title(figTitle)

    if !isempty(savefigure)
        PyPlot.savefig(savefigure)
    end
    PyPlot.display_figs()


end
