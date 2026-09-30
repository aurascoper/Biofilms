# Radial and contact follow-ups completed 2026-09-29

This record reports the completed pre-registered coefficient, horizon and
same-species-contact analyses. It does not change the biological scope of the
model. The outputs are computational results about the shipped CPM acceptance
rule.

## Fixed run identity

- Source revision: `66ffe9fca8b463dad83088dc501daf23e1d7ed3d`
- Run archive: `/home/aurascoper/Developer/Biofilms-runs/2026-09-29_66ffe9f`
- Seeds: 42--297, paired within seed
- Configuration: `N=40`, six parcels per species
- Radial shards validated: 4,352 of 4,352
- Contact shards validated: 16 of 16

Every radial shard passed the registered provenance, completeness, baseline,
shape and finite-value checks. Every contact shard also passed row-level contact
reconstruction and the archived-radius comparison.

## Registered results

The direct-radiation vector sweep did not reach the registered inward criterion
at any tested multiplier from 1 through `10^6`. At `10^6`, the
*C. neoformans* paired mean for scaled full minus null was
`-0.013697793 r/R`, with a normal 95% interval of
`[-0.0171072, -0.0102884]`.

The melanin-production vector sweep first reached the registered criterion at
`10^4`. The *C. neoformans* paired mean was `-0.081488789 r/R`, with a
normal 95% interval of `[-0.0877778, -0.0751998]`. The crossing did not remain
at `10^6`, where the mean was `-0.003125227 r/R`. Both sweeps were non-monotone
on the registered decade grid. The first crossing is not a continuous threshold
or a calibrated coefficient.

At the shipped coefficients and 400 MCS, the *C. neoformans* full-minus-null
paired mean was `+0.001496547 r/R`, with a normal 95% interval of
`[-0.0013332, +0.0043263]`. It remained inside the declared `+/-0.05 r/R`
band.

The contact prediction failed. Full minus melanin-production-off gave a mean
same-species boundary-contact difference of `-0.00130668` for
*C. neoformans* and `-0.00348245` for *C. sphaerospermum*. Their one-sided
exact sign-test p-values were `0.91270892` and `0.72389646`. The paired-radius
90% intervals remained inside `+/-0.05 r/R`, but neither contact effect met the
registered positive direction or `0.05` minimum size.

## Receipt hashes

| File | SHA-256 |
|---|---|
| `radial/SUMMARY.json` | `5aac75bd9318bd4266f2a36578289b7861d7d203ed13ff70ac98966bd21d3fdb` |
| `radial/ANALYSIS_RECEIPT.json` | `a326ff60862b665c2e5f6ddb13fff42dc5076b47602334de8c25950ca5ebebcb` |
| `radial/BASELINE_VALIDATED.json` | `30fb8289a5740338bc3ba91cdc81a6252539e3c8cfbdceb6829b4a00a4e03ddb` |
| `radial/QUEUE_PLAN.json` | `c3e858a364dd85b7300784514137f46ceab0bbce24dcafb49375fe0578cceb26` |
| `radial/QUEUE_STATUS.json` | `5b80a78b28a45e74a6bbadfb49c71e925ed9f9b51c0f95535579cf84a7af47ad` |
| `clustering/SUMMARY.txt` | `33bd603f1fb0c549f66837be04d3dd57160762c0a2660a978c9af562dd4902fa` |
| `clustering/ANALYSIS_RECEIPT.json` | `66e6b26dc1c9ffd30e66852a5b2c7c1025681aacf4bd2d80ebf6b190e1f26b21` |
| `clustering/QUEUE_PLAN.json` | `f4368ded9d53699797e8129cada6a42a70ecb8596b200cb7a4ca59222d7de630` |
| `clustering/QUEUE_STATUS.json` | `6b9ad4385f0b5536d7b2e5cdb815077288aef3d6c1a35b1da2634f532a50f18e` |

The radial summary and contact summary were regenerated from the archived
shards after completion and compared byte for byte with the recorded outputs.
