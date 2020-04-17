include("sirx.jl")
include("evolution.jl")

# Hubei
N=57.1e6
d = getParams(0.038,0.073,6.2,8.0)
s0 = initialize(N,0.9*10.0^3,2.55,d)
(s,ds) = evolve(N,s0,d,100);
plotEvolution(N,s,d,log10=true)
