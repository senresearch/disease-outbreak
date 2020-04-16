struct SIRXDynamics
    α::Float64
    β::Float64
    κ::Float64
    κ0::Float64
end

function initialize(N::Float64,C0::Float64,IXRatio::Float64,d::SIRXDynamics)
    state = zeros(nstates(d))
    state[4] = C0/N # X0, page 1, Supplement
    state[2] = IXRatio*(C0/N) # I0, page 1, supplement
    state[1] = 1.0-state[2]-state[4]
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
    S = -( d.α * st.I + d.κ0) * st.S
    I = d.α*st.I*st.S - d.β*st.I - d.κ0*st.I - d.κ*st.I
    R = d.β*st.I + d.κ0*st.S
    X = (d.κ+d.κ0) * st.I
    return [S,I,R,X]
end

function nstates(d::SIRXDynamics)
    return 4
end
function stateNames(d::SIRXDynamics)
    return ["S" "I" "R" "X"]
end

function getParams(κ::Float64,κ0::Float64,
    R0Free::Float64,TInfected::Float64)
    β = 1.0/TInfected
    α = R0Free*β
    return SIRXDynamics(α,β,κ,κ0)
end
