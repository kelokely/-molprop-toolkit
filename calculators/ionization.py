#!/usr/bin/env python3
"""Ionization / pH-dependent helper utilities.

This module supports two conceptual approaches to pH handling.

Heuristic mode is always available and produces a small set of pH-dependent
features derived from estimated pKa values using Henderson–Hasselbalch. This
mode is designed to be lightweight and dependency-free, and it does not attempt
to generate a new, chemically enumerated protonation-state SMILES.

Enumeration mode (Dimorphite-DL) is implemented separately as an optional
feature and is added in a follow-up subtask so that this module can remain
usable even in environments where dimorphite_dl is not installed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional


@dataclass
class IonizationConfig:
    ph: float = 7.4


def frac_protonated_base(pka: float, ph: float) -> float:
    """Fraction protonated for a base BH+ <-> B + H+.

    ratio BH+/B = 10^(pKa - pH)
    fraction(BH+) = 1 / (1 + 10^(pH - pKa))
    """

    return 1.0 / (1.0 + 10 ** (ph - pka))


def frac_deprotonated_acid(pka: float, ph: float) -> float:
    """Fraction deprotonated for an acid HA <-> A- + H+.

    ratio A-/HA = 10^(pH - pKa)
    fraction(A-) = 1 / (1 + 10^(pKa - pH))
    """

    return 1.0 / (1.0 + 10 ** (pka - ph))


def estimate_net_charge(basic_pka: Optional[float], acidic_pka: Optional[float], ph: float) -> float:
    """Estimate net charge using the simplest one-acid/one-base model.

    This is a coarse heuristic that is mainly useful as a ranking signal and for
    selecting among enumerated protomers.
    """

    z = 0.0
    if basic_pka is not None:
        z += frac_protonated_base(basic_pka, ph)
    if acidic_pka is not None:
        z -= frac_deprotonated_acid(acidic_pka, ph)
    return z


def classify_state(basic_pka: Optional[float], acidic_pka: Optional[float], ph: float) -> str:
    """Classify dominant ionization type at a given pH."""

    fp = frac_protonated_base(basic_pka, ph) if basic_pka is not None else 0.0
    fd = frac_deprotonated_acid(acidic_pka, ph) if acidic_pka is not None else 0.0

    if basic_pka is not None and acidic_pka is not None and fp > 0.5 and fd > 0.5:
        return "zwitterionic"

    z = fp - fd
    if z > 0.5:
        return "cationic"
    if z < -0.5:
        return "anionic"
    return "neutral"


def estimate_logd(logp: float, basic_pka: Optional[float], acidic_pka: Optional[float], ph: float) -> float:
    """Approximate logD at a given pH from logP and pKa.

    This uses the same single-site approximation commonly used in quick triage.
    It is not a replacement for measured logD or more detailed microspecies
    modeling.
    """

    if basic_pka is not None and basic_pka > ph:
        return logp - math.log10(1 + 10 ** (basic_pka - ph))
    if acidic_pka is not None and acidic_pka < ph:
        return logp - math.log10(1 + 10 ** (ph - acidic_pka))
    return logp


def ionization_features(
    logp: float,
    basic_pka: Optional[float],
    acidic_pka: Optional[float],
    ph: float = 7.4,
) -> Dict[str, object]:
    """Return a small set of pH-dependent ionization features."""

    fp = frac_protonated_base(basic_pka, ph) if basic_pka is not None else None
    fd = frac_deprotonated_acid(acidic_pka, ph) if acidic_pka is not None else None
    z = estimate_net_charge(basic_pka, acidic_pka, ph)
    state = classify_state(basic_pka, acidic_pka, ph)
    logd = estimate_logd(logp, basic_pka, acidic_pka, ph)

    out: Dict[str, object] = {
        "Ion_pH": float(ph),
        "Ion_State": state,
        "Ion_NetCharge_Est": round(z, 3),
        "Ion_LogD_pH_Est": round(logd, 3),
    }

    if fp is not None:
        out["Ion_FracProtonatedBase"] = round(fp, 3)
    if fd is not None:
        out["Ion_FracDeprotonatedAcid"] = round(fd, 3)

    # Convenience integer for selection logic; this is intentionally rough.
    out["Ion_NetCharge_Round"] = int(round(z))

    return out

