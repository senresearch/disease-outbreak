include("seir.jl")

d0 = getParams(5.0,6.0,0.5,0.1,0.1,0.15,0.06,6.0,8.0,0.4)
d1 = getParams(5.0,6.0,0.5/2,0.1/2,0.1/2,0.15,0.06,6.0,8.0,0.4)

state0 = initialize(1000.0,1.0,d0)
(states0,deltaStates0) = evolve(state0,d0,300)
plotEvolution(states0,d0)

state0 = initialize(1000.0,1.0,d0)
(states0,deltaStates0) = evolve(state0,d0,60)
(states1,deltaStates1) = evolve(states0[:,60],d1,241)
plotEvolution([states0 states1[:,2:end]],d0)

state0 = initialize(1000.0,1.0,d0)
(states0,deltaStates0) = evolve(state0,d0,30)
(states1,deltaStates1) = evolve(states0[:,30],d1,271)
plotEvolution([states0 states1[:,2:end],d0)

state1 = initialize(1000.0,1.0,d1)
(states1,deltaStates1) = evolve(state1,d1,300)
plotEvolution(states1,d1)
