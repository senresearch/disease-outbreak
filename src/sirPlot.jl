"""
sirPlot
========

Synopsis
--------
*sirPlot* is a collection of plotting function related to the different SIR
models.

Exported Functions
------------------
* ``plotEvolution``: Plots evolution of the states for the SIR model.

* ``plotEvoFitted``: Plots the Activity data array over 24 hours with two plots on
                   the same axes with different left and right scales.



See also
--------


References
----------

"""

#include("ActStatData.jl")

using PyCall, PyPlot
pygui(:tk)

include(realpath(joinpath((@__DIR__), "sirx.jl")))


#-------------------------------------------------------------------------------
"""
plotEvolution(N::Float64, states::Matrix{Float64},
              d::Dynamics;  log10 = "Null", figTitle = "",
              grid = "off")

Generate a plot of states evolved for the SIR model.

# Arguments
- N = the population size
- states = a matrix of states evoloved over time. See the function, `evolve`.
- d = a type of Dynamics of the SIR model
- log10 = a keyword argument (default is "Null") where "semilog" makes a plot
          with log scaling on the y axis, "loglog" makes a plot with log scaling
          on the x axis and y axis, and "Null"makes a plot without log scaling.
- figTitle = a keyword argument (default is empty) that add a title to the figure.
- grid = a keyword argument (default is `false`) that adds a grid on the plotting if
         it is `true`.
- fsize = a keyword argument (default is [5, 5]]) that set the size of the figure.
- savefigure = a keyword argument (default is empty "") that save the figure in
               a given path.(e.g., savefigure = "../images/Plot.png") The extension
               of the filename indicates the format that is used to save the figure
               such as .pdf, .png, .jpg,...
"""

function myPlotEvolution(N::Float64, states::Matrix{Float64},
                         d::Dynamics; log10 = "Null", figTitle = "",
                         grid = false, fsize= [5, 5], savefigure = "")

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
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
"""
plotEvoFitted(N::Float64, states::Matrix{Float64},
              d::Dynamics, actualState::Array{Float64, 1};
              log10 = "Null", figTitle = "", grid = "off")

Generate a plot of fitted state evolution for the SIR model.

# Arguments
- pop = a type SIRXPopulation that contains population size
- states = a matrix of states evoloved over time. See the function, `evolve`.
- d = a type of Dynamics of the SIR model
- actualState = actual case numbers
- log10 = a keyword argument (default is "Null") where "semilog" makes a plot
          with log scaling on the y axis, "loglog" makes a plot with log scaling
          on the x axis and y axis, and "Null"makes a plot without log scaling.
- figTitle = a keyword argument (default is empty) that add a title to the figure.
- grid = a keyword argument (default is `false`) that adds a grid on the plotting if
         it is `true`.
- fsize = a keyword argument (default is [5, 5]]) that set the size of the figure.
- savefigure = a keyword argument (default is empty "") that save the figure in
               a given path.(e.g., savefigure = "../images/Plot.png") The extension
               of the filename indicates the format that is used to save the figure
               such as .pdf, .png, .jpg,...
"""

function plotEvoFitted(N::Float64, states::Matrix{Float64},
                         d::Dynamics, actualState::Array{Float64, 1};
                         log10 = "Null", figTitle = "", grid = false,
                         fsize= [5, 5], savefigure = "")


    fig, (ax) = PyPlot.subplots(nrows=1, figsize=(fsize[1], fsize[2]))

    ###############
    #  PLOT DATA  #
    ###############





    if log10 == "semilog"
        ax.semilogy(N*states'[:, 4], label=stateNames(d)[4]*"ₜ fitted")
        ax.semilogy(actualState, linestyle ="-", color = "#A7A7A7",
                  marker="s", markersize = 15, markerfacecolor= "none",
                  label=stateNames(d)[4]*"ₜ  actual")
        ax.semilogy(N*states'[:, 2], linestyle ="--", color = "red",
                  label=stateNames(d)[2]*"ₜ  fitted")
    elseif log10 == "loglog"
        ax.loglog(N*states'[:, 4], label=stateNames(d)[4]*"ₜ  fitted")
        ax.loglog(actualState, linestyle ="-", color = "#A7A7A7",
                  marker="s", markersize = 15, markerfacecolor= "none",
                  label=stateNames(d)[4]*"ₜ  actual")
        ax.loglog(N*states'[:, 2], linestyle ="--", color = "red",
                  label=stateNames(d)[2]*"ₜ  fitted")
        ax.set_ylim(1.0e1, 1.0e5)
    elseif log10 == "null"
        ax.plot(N*states'[:, 4], label=stateNames(d)[4]*"ₜ  fitted")
        ax.plot(actualState, linestyle ="-", color = "#A7A7A7",
                  marker="s", markersize = 15, markerfacecolor= "none",
                  label=stateNames(d)[4]*"ₜ  actual")
        ax.plot(N*states'[:, 2], linestyle ="--", color = "red",
                  label=stateNames(d)[2]*"ₜ  fitted")
        ylim(-3000, maximum(tt)*1.10)
        xlim(-1, 40)
    else
        println("Incorrect log10 option, no log is applied.")
        ax.plot(N*states'[:, 4], label=stateNames(d)[4]*"ₜ  fitted")
        ax.plot(actualState, linestyle ="-", color = "#A7A7A7",
                  marker="s", markersize = 15, markerfacecolor= "none",
                  label=stateNames(d)[4]*"ₜ  actual")
        ax.plot(N*states'[:, 2], linestyle ="--", color = "red",
                  label=stateNames(d)[2]*"ₜ  fitted")
    end


    #####################
    #  AXES AND LABELS  #
    #####################

    ax.spines["top"].set_visible(false) # Hide the top edge of the axis
    ax.spines["right"].set_visible(false) # Hide the right edge of the axis
    ax.xaxis.set_ticks_position("bottom")
    # ax.spines["left"].set_position(("axes",-0.03)) # Offset the left scale from the axis
    # ax.spines["bottom"].set_position(("axes",-0.05)) # Offset the bottom scale from the axis



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
#-------------------------------------------------------------------------------
