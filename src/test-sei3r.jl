include("sei3r.jl")
include("evolution.jl")

# default parameters for RCAN model
d0 = getParams(5.0,6.0,0.5,0.1,0.1,0.15,0.06,6.0,8.0,0.4,SEI3R())
# social distancing reducing transmission probs by 1/2
d1 = getParams(5.0,6.0,0.5/2,0.1/2,0.1/2,0.15,0.06,6.0,8.0,0.4,SEI3R())

N = 1000.0
state0 = initialize(1.0/N,d0)
(states0,deltaStates0) = evolve(N,state0,d0,300)
plotEvolution(N,states0,d0)
# savefig("dist300.pdf")
state0 = initialize(1.0/N,d0)
(states0,deltaStates0) = evolve(N,state0,d0,60)
(states1,deltaStates1) = evolve(N,states0[:,60],d1,241)
plotEvolution(N,[states0 states1[:,2:end]],d0)
# savefig("dist60.pdf")
state0 = initialize(1.0/N,d0)
(states0,deltaStates0) = evolve(N,state0,d0,30)
(states1,deltaStates1) = evolve(N,states0[:,30],d1,271)
plotEvolution(N,[states0 states1[:,2:end]],d0)
# savefig("dist30.pdf")
state1 = initialize(1.0/N,d1)
(states1,deltaStates1) = evolve(N,state1,d1,300)
plotEvolution(N,states1,d1)
# savefig("dist0.pdf")
