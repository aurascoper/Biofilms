# Brief 3: evidence gates for the phase passenger

This brief asks whether the prospective phase passenger is worth exercising. It does not change the model, supply biological parameters, or turn the phase specification into a reported method.

Status: research brief, 2026-09-29. One precedent gate is partly closed, one theory calculation is reproduced below, and the biological parameter search remains open.

## 1. Liu 2017 confirms a phase-oscillator precedent, but its exact law is still unverified

Liu et al. experimentally observed synchronized growth and electrical-signalling oscillations in pairs of *Bacillus subtilis* biofilms. Their main text says that they modelled the two biofilms as “two coupled phase oscillators” and used the model to predict in-phase and antiphase regions as communication and nutrient competition changed. This is an empirical precedent for treating an entire biofilm as an oscillator.

The accessible main article does not print the phase equations. It points to its supplementary text. The supplementary PDF could not be retrieved through either the Science publisher route or the PubMed Central attachment route on 2026-09-29; the institutional repository copy contains the main article only. The evidence therefore does not yet support the narrower claim that Liu used a Kuramoto sine-coupling law. That claim remains open until the supplementary equation is retrieved, quoted and archived with its definitions.

Even if the equation is Kuramoto-type, the unit of the precedent is one whole biofilm. The proposed passenger assigns one oscillator to each CPM parcel on a moving contact graph. That is an implementer's analogue, not a reproduction of Liu's biological model.

Primary article: <https://doi.org/10.1126/science.aah4204>. Accessible manuscript record: <https://pmc.ncbi.nlm.nih.gov/articles/PMC5645014/>.

## 2. Biological searches that remain open

Each search returns a source-backed value, an explicit absence after a documented search, or an unresolved retrieval failure. A positive claim requires the quoted method or result and the source location.

1. Find measured oscillation periods for *C. neoformans*, *D. radiodurans*, *C. sphaerospermum*, *A. niger*, *S. oneidensis* and *O. intermedium* (now *Brucella intermedia* in current taxonomy). The existing *B. subtilis* period does not transfer to them.
2. Search for an ionising-radiation-induced phase shift in a microbial oscillator. Record radiation modality, dose, dose rate, oscillator, phase estimator and control. Evidence from mammalian circadian systems does not fill this gate.
3. Trace the 0.36-s-per-MCS convention to its first source. Until its physical calibration is shown, treat it as inherited rather than measured.
4. Repeat the census of cell-based models with a per-cell or per-agent oscillator. The earlier two-precedent count is provisional until the search terms, databases and cut-off date are recorded.
5. Audit the assumptions behind the all-to-all critical-coupling result, the finite-sample reference and the claimed absence of synchronization below five spatial dimensions. A result for an infinite static lattice does not automatically transfer to a finite moving graph of tens of parcels.
6. Retrieve Liu's supplementary equation. A second documented retrieval failure is an admissible result; an unsupported confirmation is not.

## 3. The normal-frequency reference reproduces the expected critical coupling

For the infinite-population, all-to-all Kuramoto model with a symmetric unimodal frequency density `g`, the nonzero stationary order parameter satisfies

`1 = K integral[-pi/2, pi/2] cos(theta)^2 g(K r sin(theta)) dtheta`.

Taking `r` to zero gives `1 = K_c g(0) pi/2`, hence `K_c = 2/[pi g(0)]`. For normally distributed intrinsic frequencies with standard deviation `sigma`, `g(0) = 1/[sigma sqrt(2 pi)]`, so

`K_c = sigma sqrt(8/pi) = 1.595769121605731 sigma`.

Direct quadrature of the same self-consistency equation with `sigma = 1` gives:

| `K/K_c` | stationary `r` |
|---:|---:|
| 0.5 | 0 |
| 0.9 | 0 |
| 1.0 | 0 at the critical limit |
| 1.1 | 0.509949445943 |
| 2.0 | 0.936985007104 |

These values describe the infinite all-to-all continuum. They do not predict onset on the finite, degree-normalised, moving contact graph in `phase_locking_spec.md`. The implementation control must use the exact finite frequency sample and must report both the order parameter and effective-frequency spread.

## 4. Decision rule

The layer is worth running first as a synthetic implementation control. It is ready for a biological interpretation only after at least one relevant oscillator period and a defensible graph and coupling law are found. Radiation coupling remains `none` until a microbial radiation-to-phase measurement supplies a sign, exposure variable and rate law.

No repository source or unpublished manuscript is needed for these literature searches. An external research session should receive only this question list and public citations.
