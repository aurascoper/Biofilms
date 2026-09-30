module SyntheticResponse

export ResponseConfig

const SPECIES = (:CN, :DR, :CS, :BS, :AN, :SO, :OI)

"""
Validated inputs for a proposed synthetic dose-response experiment.
All coefficients are explicit; this object has no stepper or automatic effect
on a CPM state. The immutable tuple prevents later edits to an input vector
from bypassing validation. Units are Gy⁻¹, seconds/MCS and Gy/second.
"""
struct ResponseConfig
    alpha_Gy_inv::NTuple{7,Float64}
    seconds_per_mcs::Float64
    axis_rate_Gy_s::Float64

    function ResponseConfig(; alpha_Gy_inv, seconds_per_mcs, axis_rate_Gy_s)
        (alpha_Gy_inv isa AbstractVector || alpha_Gy_inv isa Tuple) &&
            length(alpha_Gy_inv) == length(SPECIES) ||
            throw(ArgumentError("alpha_Gy_inv must have one entry per species in CN,DR,CS,BS,AN,SO,OI order"))
        all(x -> x isa Real && isfinite(x) && x >= 0, alpha_Gy_inv) ||
            throw(ArgumentError("alpha_Gy_inv must be finite and nonnegative"))
        alpha = Tuple(Float64.(alpha_Gy_inv))
        all(isfinite, alpha) || throw(ArgumentError("alpha_Gy_inv must fit in Float64"))
        seconds_per_mcs isa Real && isfinite(seconds_per_mcs) && seconds_per_mcs > 0 ||
            throw(ArgumentError("seconds_per_mcs must be finite and strictly positive"))
        dt = Float64(seconds_per_mcs)
        isfinite(dt) && dt > 0 || throw(ArgumentError("seconds_per_mcs must fit in positive Float64"))
        axis_rate_Gy_s isa Real && isfinite(axis_rate_Gy_s) && axis_rate_Gy_s >= 0 ||
            throw(ArgumentError("axis_rate_Gy_s must be finite and nonnegative"))
        axis = Float64(axis_rate_Gy_s)
        isfinite(axis) || throw(ArgumentError("axis_rate_Gy_s must fit in Float64"))
        new(alpha, dt, axis)
    end
end

end
