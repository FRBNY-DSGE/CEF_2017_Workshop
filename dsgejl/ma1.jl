using DataStructures

"""
```
MA1{T} <: AbstractModel{T}
```

Implements the following MA(1) model:

```
x_t = μ + u_t + α*u_{t-1}
u_t ∼ N(0, σ^2)
```
"""
type MA1{T} <: AbstractModel{T}
    parameters::ParameterVector{T}                         # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                       # model steady-state values
    keys::OrderedDict{Symbol,Int}                          # human-readable names for all the model
                                                           # parameters and steady-states

    endogenous_states::OrderedDict{Symbol,Int}             # these fields used to create matrices in the
    exogenous_shocks::OrderedDict{Symbol,Int}              # measurement and equilibrium condition equations.
    expected_shocks::OrderedDict{Symbol,Int}               #
    equilibrium_conditions::OrderedDict{Symbol,Int}        #
    endogenous_states_augmented::OrderedDict{Symbol,Int}   #
    observables::OrderedDict{Symbol,Int}                   #

    spec::String                                           # Model specification number (eg "m990")
    subspec::String                                        # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                         # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}                    # Settings/flags for testing mode
    rng::MersenneTwister                                   # Random number generator
    testing::Bool                                          # Whether we are in testing mode or not

    observable_mappings::OrderedDict{Symbol, Observable}
end

DSGE.description(m::MA1) = "Moving average model of order 1: MA1, $(m.subspec)"

function MA1(subspec::String="ss0"; testing = false)
    # Model-specific specifications
    spec               = "ma1"
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister()

    # initialize empty model
    m = MA1{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(), OrderedDict{Symbol,Int}(),

            # model indices
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}())

    # Set settings
    settings_ma1!(m)
    default_test_settings!(m)

    # Initialize parameters
    init_parameters!(m)

    init_model_indices!(m)
    steadystate!(m)

    return m
end

"""
```
init_parameters!(m::MA1)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::MA1)
    m <= parameter(:μ, 0.5, (-1e5, 1e5), (1e-5, 5.0), DSGE.Exponential(),
                   DSGE.Normal(0.0, 1.0), fixed = false,
                   description = "μ: constant coefficient",
                   tex_label = "\\mu")
    m <= parameter(:β, 0.5, (-1e5, 1e5), (1e-5, 5.0), DSGE.Exponential(),
                   DSGE.Normal(0.0, 1.0), fixed = false,
                   description = "β: coefficient on u_{t-1}",
                   tex_label = "\\beta")
    m <= parameter(:σ, 0.1, (1e-5, 5.0), (1e-5, 5.0), DSGE.Exponential(),
                   DSGE.RootInverseGamma(5.0, 0.5), fixed = false,
                   description = "σ: standard deviation of u_t",
                   tex_label = "\\sigma")
end

"""
```
init_model_indices!(m::MA1)
```

Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::MA1)
    # Endogenous states
    endogenous_states = [:u_t, :u_t1]

    # Exogenous shocks
    exogenous_shocks = [:u_t]

    # Expectations shocks
    expected_shocks = Symbol[]

    # Equilibrium conditions
    equilibrium_conditions = [:eq_u_t, :eq_u_t1]

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = []

    # Measurement equation observables
    observables = [:x_t]

    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(exogenous_shocks);            m.exogenous_shocks[k]            = i end
    for (i,k) in enumerate(expected_shocks);             m.expected_shocks[k]             = i end
    for (i,k) in enumerate(equilibrium_conditions);      m.equilibrium_conditions[k]      = i end
    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(observables);                 m.observables[k]                 = i end
end

"""
```
steadystate!(m::MA1)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function DSGE.steadystate!(m::MA1)
    return m
end

function settings_ma1!(m::MA1)
    default_settings!(m)

    m <= Setting(:use_population_forecast, false, "Whether to use population forecasts as data")
end

function DSGE.eqcond(m::MA1)
    endo = m.endogenous_states
    exo  = m.exogenous_shocks
    eq   = m.equilibrium_conditions

    Γ0 = zeros(n_states(m), n_states(m))
    Γ1 = zeros(n_states(m), n_states(m))
    C  = zeros(n_states(m))
    Ψ  = zeros(n_states(m), n_shocks_exogenous(m))
    Π  = zeros(n_states(m), n_shocks_expectational(m))

    # Γ0*s_t = Γ1*s_{t-1} + Ψ*ϵ_t + Π*η_t + C
    # s_t = [u_t, u_{t-1}]'
    # ϵ_t = [u_t]'

    # Row 1: u_t = 0*u_{t-1} + 0*u_{t-2} + u_t
    Γ0[eq[:eq_u_t], endo[:u_t]] = 1
     Ψ[eq[:eq_u_t],  exo[:u_t]] = 1

    # Row 2: u_{t-1} = u_{t-1} + 0 u_{t-2} + 0 u_t
    Γ0[eq[:eq_u_t1], endo[:u_t1]] = 1
    Γ1[eq[:eq_u_t1], endo[:u_t]] = 1

    return Γ0, Γ1, C, Ψ, Π
end

function DSGE.measurement{T<:AbstractFloat}(m::MA1{T}, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T};
                                            shocks::Bool = false)
    endo = m.endogenous_states # OrderedDict{Symbol, Int} mapping state names (e.g. `:u_t`) to indices
    exo  = m.exogenous_shocks  # ... mapping shock names to indices
    obs  = m.observables       # ... mapping observable names to indices

    ZZ = zeros(n_observables(m), n_states(m))
    DD = zeros(n_observables(m))
    MM = zeros(n_observables(m), n_shocks_exogenous(m))
    EE = zeros(n_observables(m), n_observables(m))
    QQ = zeros(n_shocks_exogenous(m), n_shocks_exogenous(m))

    # y_t = Z*s_t + D
    # y_t = [x_t]'
    # s_t = [u_t, u_{t-1}]'

    # TODO: fill in entries of ZZ matrix
    # x_t = μ + u_t + β*u_{t-1}
    # DD[obs[:x_t]] = m[:μ]
    # ZZ[obs[:x_t], ...] = ...

    DD[obs[:x_t]] = m[:μ]
    ZZ[obs[:x_t], endo[:u_t]]  = 1
    ZZ[obs[:x_t], endo[:u_t1]] = m[:β]

    # TODO: fill in entries of QQ matrix
    # QQ[exo[:u_t], exo[:u_t]] = ...

    QQ[exo[:u_t], exo[:u_t]] = m[:σ]^2

    HH    = EE + MM*QQ*MM'
    VV    = QQ*MM'
    VVall = [[RRR*QQ*RRR' RRR*VV];
             [VV'*RRR'    HH]]

    return Measurement(ZZ, DD, QQ, EE, MM, VVall)
end

# This is an MA(1)-specific version of `solve` which we can use because MA(1)
# models are trivially rational expectations models
function DSGE.solve(m::MA1)
    Γ0, Γ1, C, Ψ, Π = eqcond(m)

    # Check that Γ0 is the identity
    @assert Γ0 == eye(n_states(m))

    # Check that Π is empty
    @assert isempty(Π)

    T = Γ1
    R = Ψ

    return T, R, C
end