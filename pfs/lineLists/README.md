# Line lists

Reference lines for the PFS calibration lamps and the night sky, read by
`pfs.drp.stella.ReadLineListTask`. The column layout and the meaning of the
status and source bits are given in each file's header (written by
`ReferenceLineSet.writeLineList`).

## One file per lamp

A line's intensity, its flags and even which lines appear depend on the lamp:
its fill gas and pressure, its impurities, and how warm it is. So the lists
are per **lamp**, not per element. `ReadLineListConfig.lightSourceMap` maps
each lamp to its file. The lines are then filtered to the elements that the
lamp contains (`lsst.obs.pfs.utils.getLampElements`).

| lamp (header) | fibres | file | elements kept |
|---|---|---|---|
| `Ar`, `Kr`, `Ne`, `Xe` | PFI (science) | `Ar.txt`, `Kr.txt`, `Ne.txt`, `Xe.txt` | the element |
| `HgCd` | PFI | `HgCd.txt` | Hg, Cd, Ar (the Ar buffer gas) |
| `HgAr` | PFI (used 2020–2024) | `HgAr.txt` | Hg, Ar |
| `<lamp>_eng` | engineering (IIS) | `<lamp>_IIS.txt` | the element(s); `He_eng` also O |

An IIS lamp (`<lamp>_eng`) uses its own `<lamp>_IIS.txt` if the map has an
entry for it. Otherwise it falls back to the PFI lamp's file. When a PFI lamp
and an IIS lamp of the same element are both lit, the PFI copy of each
shared line is kept.

`contamination.txt` (O I, via the `OI` light source), `skyLines*.txt` and
`galaxy.txt` aren't lamps.

## Intensities

Intensities flagged `EMPIRICAL_INTENSITY` (source bit 64) were measured from
PFS arc spectra. They are in counts per second, per pfsArm exposure time, and
quartz-normalised. Unflagged intensities are catalogue values, in arbitrary
units.

| file | visit | reduced in | notes |
|---|---|---|---|
| `Ar.txt` | 148314 | `u/rhl/shaka-arc` | PIPE2D-1890 |
| `Kr.txt` | 148320 | `u/rhl/shaka-arc` | PIPE2D-1890 |
| `Ne.txt` | 148326 | `u/rhl/shaka-arc` | PIPE2D-1890 |
| `Xe.txt` | 148317 | `u/rhl/shaka-arc` | PIPE2D-1890 |
| `HgCd.txt` | 148323 | `u/rhl/shaka-arc` | PIPE2D-1958; last exposure of a 3-exposure sequence |
| `He_IIS.txt` | 148332 | `u/rhl/PIPE2D-1958-He` | PIPE2D-1958; engineering fibres |
| `HgAr_IIS.txt` | 148142 | `u/rhl/iis-arcs` | PIPE2D-1958/1959; engineering fibres |
| `Ar_IIS.txt`, `Kr_IIS.txt`, `Ne_IIS.txt` | 148141, 149798, 149797 | `u/rhl/iis-arcs` | PIPE2D-1959; engineering fibres |
| `HgAr.txt`, `H_IIS.txt` | — | — | catalogue values |

How the intensities were measured:

- **Method.** All a lamp's lines are fitted at once, with known wavelengths,
  the pipeline LSF and a spline continuum, on the clipped mean of the
  quartz-normalised spectra. Faint and undetected lines are kept, with their
  errors, so they can be cut on S/N downstream.
- **Exposure time.** The fluxes are divided by the mean of the pfsArms'
  exposure times. `pfsMerged.identity.expTime` is their sum (PIPE2D-1957).
- **PFI and IIS scales differ.** The IIS (engineering-fibre) intensities are on
  their own scale and are not comparable with the PFI lamps'.
- **Bright IIS lines use integrated fluxes.** The PFI-calibrated LSF doesn't
  describe the 60 µm IIS fibres (the PFI fibres are about 120 µm), so IIS line
  groups with S/N > 30 use the flux integrated over ±0.6 nm. Blended lines are
  split by their fitted ratios.
- **Unmeasured lines.** In `HgCd.txt` and the IIS lists, lines outside the PFS
  range or otherwise unmeasured get their previous intensity, scaled by the
  median measured/previous ratio, and have no error. In Ar, Kr, Ne and Xe.txt
  they keep their catalogue values.
- **Lamp warm-up.** Lamps change as they warm up: Cd brightens about 1.6×
  relative to Hg over the first few minutes, and He is about 8% fainter in its
  first exposure. This hasn't been characterised yet.
- **Future refit.** The intensities will be refit with the PFS LSF model once
  it includes the IIS/PFI fibre sizes and the halo around each line image.

## Wavelengths and transitions

Wavelengths are in vacuum. Transitions are from NIST ASD (source bit `NIST`),
written as `lower->upper` levels in the form `<valence>_<term>_<J>`. Lines that
NIST lists as unresolved blends of several components join the tokens with
`;`. Components closer than 0.06 nm are merged into one line, at the
g·A-weighted mean Ritz wavelength, with the J values joined by `/` (e.g.
`2p_3P*_0/1/2->3d_3D_1/2/3`). `UNKNOWN` means NIST doesn't classify the line.

Some lines aren't in NIST's line list. **Cd I 1039.748** (5s5p ¹P₁ – 5s6s ¹S₀)
uses its Ritz wavelength from NIST's energy levels. Three Cd I lines ending on
5f ³F (1126.933, 1163.489, 1165.900) use PFS-measured wavelengths (±0.03 nm;
source `MANUAL|NIST`).

## BROAD

`BROAD` (status bit 16) marks lines that are pressure-broadened in **these**
lamps. A line is broadened in this way if its lower level is a J = 1
resonance level, i.e. one that decays directly to the ground state. Collisions
with the lamp's ground-state atoms broaden that level (resonance broadening),
by an amount proportional to the density times the resonance line's f·λ. So
**BROAD depends on the lamp's fill pressure as well as on the line**, and the
flags apply to the lamps that were measured. The PFI lamps' fill pressures are
Ne 3, Ar 45, Kr 60 and Xe 20 Torr.

Flags are set by family: every line ending on a broadened level is flagged,
in every arm. The excess Lorentzian width (FWHM) was measured in the r arm,
relative to the same lamp's other lines of similar brightness:

| lamp | 1s₂ (`ns_2[1/2]*_1`) | 1s₄ (`ns_2[3/2]*_1`) | flagged |
|---|---|---|---|
| Ar | +5.1 GHz (+4.2 in n) | +0.1 GHz | 1s₂ |
| Kr | +6.2 GHz | +8.1 GHz | 1s₂, 1s₄ |
| Xe | +2.9 ± 2.4 GHz | +3.3 GHz | 1s₂, 1s₄ |
| Ne | +1.4 GHz (2σ) | −0.6 GHz | none |

- **Thresholds.** A family is flagged if its excess is ≥ 3σ and ≥ 2 GHz.
- **Xe 1s₂** is flagged because its resonance line's f·λ predicts about 0.65×
  the 1s₄ width, consistent with the measurement.
- **Ne** isn't broadened at 3 Torr (about 0.1 GHz predicted). **Ar 1s₄**'s
  resonance line has about 4× less oscillator strength than 1s₂'s.
- **The b arm can't measure these widths** (about 100 GHz per pixel), so its
  flags follow from the family (PIPE2D-1945).

**The broadened families are also shifted to the blue.** In the r arm (all
four spectrographs) their centroids are offset from the same lamp's other
lines by:

| family | shift |
|---|---|
| Ar 1s₂ | −0.017 ± 0.004 px (−1.4 pm) |
| Kr 1s₂ | −0.015 ± 0.003 px |
| Kr 1s₄ | −0.018 ± 0.002 px (−1.6 pm) |
| Xe 1s₄ | −0.012 ± 0.006 px |

The unbroadened families (Ar 1s₄, Ne) show ≤ 0.004 px. The shift is about
−0.1 times the broadening, and is not in NIST's low-pressure wavelengths. The
lamps' overall offsets are small: relative to Ne, Ar is +0.008 ± 0.003 px and
Xe is 0.000 ± 0.003 px. Lines flagged BROAD should therefore not be used for
wavelength calibration at the 0.01 px level.

The IIS lamps were measured the same way, using all 64 engineering fibres
(arcs 149797, 149798 and 148141), r arm:

| lamp | 1s₂ | 1s₄ | flagged |
|---|---|---|---|
| Ar_IIS | +6.0 GHz (5σ) | −0.1 GHz | 1s₂ |
| Kr_IIS | +1.3 ± 1.7 GHz | +2.9 ± 5 GHz | none |
| Ne_IIS | +1.3 GHz (1σ) | −0.2 GHz | none |

The IIS Kr lamp's 1s₂ width is ~3σ below the PFI Kr lamp's, suggesting a
lower fill pressure. HgAr_IIS (a new lamp) and HgCd's Ar buffer-gas lines are
not flagged. `ReadLineListTask` keeps BROAD lines. They are left out of the detectorMap
fits (`FitDetectorMapTask` and `AdjustDetectorMapTask` ignore `lineFlags =
["BAD"]`, and BAD includes BROAD), the LSF fit (`FitLsfTask`,
`excludeCatalogueBroad`) and the aperture-correction fit.
