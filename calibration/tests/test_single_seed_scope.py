"""The melanin numbers are one seed's, and a published sentence that states them must say so.

M = 1.44 is C. sphaerospermum's field value in the seed-42 run, and 1.155, 15.5% and the 9.6
Delta H ratio are computed from it. Over seeds 42-297 at the same configuration that value is the
31st lowest of 256, and the order CS > CN > AN it was quoted with fails in 19 (PP-65-06). The
manuscript said "the ordering follows the hand-specified production scales" all the same, and
set 15.5% against the one-role radiation value while PP-62-13 recorded the reach as stated.
test_withdrawn_comparisons.py could not see that sentence: it matches the four-orders wording,
and the sentence had none.

WHAT BOUNDS IT. Three documents; four numbers, one ordering phrase and one comparison, each
judged inside a window. A seed named near an unrelated number satisfies it, and the same claim
in other words is invisible to it, as it is to the scan beside it.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from test_withdrawn_comparisons import flatten

REPO = Path(__file__).resolve().parents[2]
DOCS = ("preprint/modeling_radioresistance_and_radiotropic_fitness.tex", "README.md",
        "docs/calibration/integration_contract.md")

SEED_VALUE = re.compile(r"(?<![\d.,])(?:1\.44|1\.155|15\.5)(?![\d,])|(?<![\d.])9\.6 in")
# "seeds 42-297" is a range, not a scope, and must not excuse a seed-42 number.
SEED = re.compile(r"seed[~\s-]*42")
ORDERING = re.compile(r"order(?:ing)?\b[^.]{0,60}?\b(?:follows|is set by)\b"
                      r"(?=[^.]{0,60}(?:production scales|α_m|alpha_m))")
ONE_ROLE = re.compile(r"one part in (?:\$?10\^\{?5|10⁵|1e5)")
REACH = re.compile(r"\breach\b")          # not "reaches the dynamics"
CORRECTED = re.compile(r"corrected")

WINDOW = 160


def defects(text: str) -> list[str]:
    flat = flatten(text)

    def near(m):
        return flat[max(0, m.start() - WINDOW):m.end() + WINDOW]

    out = [m.group() for m in SEED_VALUE.finditer(flat) if not SEED.search(near(m))]
    out += ["ordering" for m in ORDERING.finditer(flat) if not CORRECTED.search(near(m))]
    out += ["one-role" for m in ONE_ROLE.finditer(flat)
            if "melanin" in near(m) and not REACH.search(near(m))]
    return out


def test_the_published_documents_scope_every_single_seed_claim():
    """THE ONE THAT MATTERS."""
    found = {doc: defects((REPO / doc).read_text(encoding="utf-8")) for doc in DOCS}
    assert not {doc: d for doc, d in found.items() if d}, found


# THE CONTROLS ARE THE SENTENCES THAT SHIPPED, as they stood at e3466f4.
ABSTRACT_WAS = r"""Auditing the CPM's radiation-linked terms against their shipped parameterisation gives two
measured results. The melanin term changes Metropolis acceptance by approximately 15.5\%, whereas
the direct per-species radiation term changes it by approximately one part in $10^5$ for the two
negatively signed species; only the product of melanin production and melanin coupling reaches
the dynamics, so the two are not separately identifiable."""

BENEATH_TABLE_3_WAS = r"""$\Delta H_{\mathrm{mel}}$ at $M=1.44$ & $-0.720$ & $\mathbf{1.155}$ \\
\bottomrule
\end{tabular}
\end{table}

For the two negatively signed species, the direct radiation term changes acceptance by about one
part in $10^5$, whereas the melanin term changes it by 15.5\%. Radiation can therefore influence"""

SECTION_6_4_WAS = r"""The dimensionless melanin field increases over the 100-MCS run for all three producer classes
(Figure~\ref{fig:melanin}). \textit{C.~sphaerospermum} reaches the largest final field value,
$M=1.44$, followed by \textit{C.~neoformans} and \textit{A.~niger}. The ordering follows the
hand-specified production scales; these field values have no calibration to a physical melanin mass
or concentration."""

README_FIG2_WAS = """*C. sphaerospermum* reaches 1.44. The ordering among the three producers is set by the input `α_M`
(0.14 versus 0.10 and 0.065), and the linear rise is structural"""


@pytest.mark.parametrize("text, expected", [
    (ABSTRACT_WAS, ["15.5", "one-role"]),
    (BENEATH_TABLE_3_WAS, ["1.44", "1.155", "15.5", "one-role"]),
    (SECTION_6_4_WAS, ["1.44", "ordering"]),
    (README_FIG2_WAS, ["1.44", "ordering"]),
], ids=["abstract", "beneath-table-3", "section-6.4", "readme-fig-2"])
def test_each_defect_is_found_in_the_sentence_that_shipped_it(text, expected):
    assert defects(text) == expected


def test_a_seed_range_does_not_scope_a_seed_42_number():
    assert defects("the bias spans 1.134-1.292 over seeds 42-297, against 1.155") == ["1.155"]
    assert defects("at the seed-42 field value the bias is 1.155") == []
