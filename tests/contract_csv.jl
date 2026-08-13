# Byte-level regression contract: the serial model's CSV output must match the
# golden fixture generated at the pre-coupling baseline (ec8929d). Any change —
# genealogy fields, dose fields, refactors — must keep this byte-identical,
# proving the legacy stochastic path (RNG stream included) is untouched.
#
# Runs validate_serial.jl as a subprocess so the real contract path is tested,
# include_string split-marker and all.

fixture_path = joinpath(@__DIR__, "fixtures", "serial_seed42.csv")
expected = filter(l -> startswith(l, "CSV"), readlines(fixture_path))

out = read(`$(Base.julia_cmd()) --project=$REPO $(joinpath(REPO, "validate_serial.jl")) 42`, String)
got = filter(l -> startswith(l, "CSV"), split(out, '\n'))

@test !isempty(expected)
@test got == expected
