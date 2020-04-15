# Notes on SIR model

## SIR model dynamics

Let S be number susceptible, I be the number infected, and R being the number removed (death or recovery).  It is assumed that S+I+R=N, the total number of individuals in the population.

∂S = -αSI

∂I = αSI -βI

∂R =  βI

## SIR-X model dynamics

∂S = -αSI -κ₀S

∂I = αSI -βI -κ₀I -κI

∂R =  βI + κ₀S

∂X = (κ+κ₀)I
