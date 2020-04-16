include("sirx.jl")
include("evolution.jl")

# Hubei
N=57.1e6
d = getParams(0.038,0.073,6.2,8.0)
s0 = initialize(N,10.0^3,0.001,d)
(s,ds) = evolve(N,s0,d,100);
plotEvolution(N,s,log10=true)
