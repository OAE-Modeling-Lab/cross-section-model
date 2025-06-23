"""
Defining a TracerLinearEquationOfState
this extends the LinearEquationOfState to include a term associated with an arbitrary tracer
"""

import Oceananigans.BuoyancyFormulations: with_float_type
import SeawaterPolynomials: ρ′

using Oceananigans.Utils: prettysummary
using Oceananigans.BuoyancyFormulations: AbstractEquationOfState


"""
    TracerLinearEquationOfState{FT} <: AbstractEquationOfState

Linear equation of state for seawater with tracer.
"""
struct TracerLinearEquationOfState{FT} <: AbstractEquationOfState
    thermal_expansion :: FT
    haline_contraction :: FT
    tracer_contraction :: FT
    reference_density :: FT
end

Base.summary(eos::TracerLinearEquationOfState) =
    string("TracerLinearEquationOfState(thermal_expansion=", prettysummary(eos.thermal_expansion),
                               ", haline_contraction=", prettysummary(eos.haline_contraction),
                               ", tracer_contraction=", prettysummary(eos.tracer_contraction),
                               ", reference_density=", prettysummary(eos.reference_density), ")")

Base.show(io::IO, eos::TracerLinearEquationOfState) = print(io, summary(eos))

with_float_type(FT::DataType, eos::TracerLinearEquationOfState) = TracerLinearEquationOfState(convert(FT, eos.thermal_expansion),
                                                                                  convert(FT, eos.haline_contraction),
                                                                                  convert(FT, eos.tracer_contraction),
                                                                                  convert(FT, eos.reference_density))

"""
    TracerLinearEquationOfState([FT=Float64;] thermal_expansion=1.67e-4, haline_contraction=7.80e-4)

Return `LinearEquationOfState` for `SeawaterBuoyancy` with
`thermal_expansion` coefficient and `haline_contraction` coefficient.
The buoyancy perturbation ``b`` for `LinearEquationOfState` is

```math
    b = g (α T - β S - γ C),
```

where ``g`` is gravitational acceleration, ``α`` is `thermal_expansion`, ``β`` is
`haline_contraction`, ``γ`` is the `tracer_contraction`,  ``T`` is temperature, and ``S`` is practical salinity units.

Default constants in units inverse Kelvin and practical salinity units
for `thermal_expansion` and `haline_contraction`, respectively,
are taken from Table 1.2 (page 33) of Vallis, "Atmospheric and Oceanic Fluid
Dynamics: Fundamentals and Large-Scale Circulation" (2nd ed, 2017).
"""
TracerLinearEquationOfState(FT=Oceananigans.defaults.FloatType; thermal_expansion=1.67e-4, haline_contraction=7.80e-4, tracer_contraction=0.0, reference_density=1029.0) =
    TracerLinearEquationOfState{FT}(thermal_expansion, haline_contraction, tracer_contraction, reference_density)

#####
##### Thermal expansion and haline contraction coefficients
#####

@inline thermal_expansion(Θ, sᴬ, D, eos::TracerLinearEquationOfState) = eos.thermal_expansion
@inline haline_contraction(Θ, sᴬ, D, eos::TracerLinearEquationOfState) = eos.haline_contraction

# Shortcuts
@inline  thermal_expansionᶜᶜᶜ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.thermal_expansion
@inline  thermal_expansionᶠᶜᶜ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.thermal_expansion
@inline  thermal_expansionᶜᶠᶜ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.thermal_expansion
@inline  thermal_expansionᶜᶜᶠ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.thermal_expansion

@inline haline_contractionᶜᶜᶜ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.haline_contraction
@inline haline_contractionᶠᶜᶜ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.haline_contraction
@inline haline_contractionᶜᶠᶜ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.haline_contraction
@inline haline_contractionᶜᶜᶠ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.haline_contraction

@inline tracer_contractionᶜᶜᶜ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.tracer_contraction
@inline tracer_contractionᶠᶜᶜ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.tracer_contraction
@inline tracer_contractionᶜᶠᶜ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.tracer_contraction
@inline tracer_contractionᶜᶜᶠ(i, j, k, grid, eos::TracerLinearEquationOfState, C) = eos.tracer_contraction
#####
##### Convinient aliases to dispatch on
#####

const TracerLinearSeawaterBuoyancy = SeawaterBuoyancy{FT, <:TracerLinearEquationOfState} where FT
#const TracerLinearTemperatureSeawaterBuoyancy = SeawaterBuoyancy{FT, <:TracerLinearEquationOfState, <:Nothing, <:Number} where FT
#const TracerLinearSalinitySeawaterBuoyancy = SeawaterBuoyancy{FT, <:TracerLinearEquationOfState, <:Number, <:Number} where FT

#####
##### buoyancy perturbation
#####

@inline buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, b::TracerLinearSeawaterBuoyancy, C) =
    @inbounds b.gravitational_acceleration * (b.equation_of_state.thermal_expansion  * C.T[i, j, k] -
                                              b.equation_of_state.haline_contraction * C.S[i, j, k] - 
                                              b.equation_of_state.tracer_contraction   * C.c[i, j, k])



@inline ρ′(T, S, C, eos::TracerLinearEquationOfState) = 
    eos.thermal_expansion * T - eos.haline_contraction * S - eos.tracer_contraction * C

