"""Saturation-detection utility — planned post-MVP.

The core ``fit`` path accepts a pixel-level ``validMask`` on each ``Ramp``.
This module will eventually provide a default saturation detector that
produces such a mask from a ramp's raw deltas / cumulative signal. For the
MVP, it intentionally contains no implementation.

See ``docs/superpowers/specs/2026-04-16-relin-package-design.md``
section 1 ("Planned (post-MVP) extensions") and section 11.
"""

from __future__ import annotations

# Intentionally empty.
