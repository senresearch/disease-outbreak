struct SIRXDynamics
    α::Float64
    β::Float64
    κ::Float64
    κ0::Float64
end

function initialize(popSize::Float64,
    IXRatio::Float64,d::SIRXDynamics)
    state = zeros(nstates(d))
    state[2] = nInfected
    state[1] = popSize-nInfected
    return state
end

# change in a day
"""
change(s::Vector{Float64},d::SIRXDynamics)

Retutn the change in the state of the population in a day

s = state vector (S,I,R)
d = SIRX dynamics parameters
"""
function change(s::Vector{Float64},d::SIRXDynamics)
    # this is done purely for readability of the formula
    st = NamedTuple{(:S,:I,:R,:X)}(s)
    N = sum(s) # population size
    S = -(( d.α * st.I + d.κ0)/N) *st.S
    I = d.α*st.I*st.S/N - d.β*st.I - d.κ0*st.I - d.κ*st.I
    R = d.β*st.I + d.κ0*st.S/N
    X = (d.κ+d.κ0)*st.I
    return [S,E,I,R,X]
end

function nstates(d::SIRXDynamics)
    return 4
end
function stateNames(d::SIRXDynamics)
    return ["S" "I" "R" "X"]
end

function getParams(κ::Float64,κ₀::Float64,
    R₀::Float64,TInfected::Float64)
    Q = (κ+κ₀)*TInfected
    β = (κ+κ₀)/Q - (κ+κ₀)
    α = R₀*(β+κ+κ₀)
    return SIRXDynamics(α,β,κ,κ₀)
end
