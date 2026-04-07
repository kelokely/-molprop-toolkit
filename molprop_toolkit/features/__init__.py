"""Feature generation subpackage.

This namespace contains reusable, dependency-light building blocks for ML-oriented exports.
RDKit is imported lazily inside functions so the broader toolkit can be installed without RDKit
(for CSV-only analysis/reporting workflows).
"""

from __future__ import annotations
