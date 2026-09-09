# Rosette phantom pipeline — fixed test build

Self-contained copy of the rosette MRSI pipeline with all fixes applied.
Run from this folder (or add it to the MATLAB path first — its files must resolve
before the older copies elsewhere on the path).

## Files
- `run_MRSI_Rosette_40x40.m` — the pipeline, all fixes below applied.
- `remove_receiver_dc.m` — receiver DC-offset removal (FOV-center dropout fix).
- `apply_water_scaling_all.m` — 3-compartment ground truth + optional relaxation.
- `plot_coord.m` — plot an LCModel `.coord` (data/fit/residual) to inspect a voxel.

## How to run
```matlab
cd 'C:\Users\divya\Downloads\fida codes\fida27May\fida-main\test'
run_MRSI_Rosette_40x40
```

## Fixes applied (what changed and why)

1. **`smoothFwhm = 0`** (was 20). 20 mm Gaussian blends the Small(70 mL) /
   Large(650 mL) / Outer(3800 mL) compartments across their borders → wrong
   per-compartment mM. Must be 0 for quantification.

2. **Spectral zero-fill DISABLED** (`doZeroFill=false`). `op_CSIspecZeroFill`
   (576→4096) pads a non-decayed FID → a step at pt 576 → uniform sinc **ripple**
   ("squiggly lines") locked to 1/(576·dt), worst where SNR is low. Proven in
   diagnose: OFF → S/N up, residual flat, LCModel fits fine at native 576.

3. **Receiver DC removal — step 1b** (`removeDC=true`). A constant per-coil ADC
   bias is uniform across k; the rosette k=0 revisit concentrates it into the
   **FOV-center voxel** → center dropout + inverted water + LCModel CRLB 999.
   `remove_receiver_dc` estimates it from high-|k| samples (object ≈ 0 there) and
   subtracts. It PRINTS `|DC|/|k0 signal|` per coil:
   - **several %+** → DC was real, now removed, center should fill in.
   - **< ~2%** → DC innocent → center dropout is coil sensitivity / B1 (large
     phantom, surface array) → NOT a recon fix; set `removeDC=false`.

4. **TE-vs-basis check.** Basis is `...TE15...` but LCModel ECHOT is auto-derived
   from the data. If they differ, a warning fires → rebuild the basis at the real
   data TE (coupled peaks like Lac mis-fit otherwise).

5. **`apply_water_scaling_all`** — 3-compartment ground truth (from the
   final_phantom sheet, mM) baked in as `GT`; `doRelax=false` (raw LCModel mM;
   undoped/unmeasured → no guessed bias). Metab relaxation table = Act/Cr/Cho/Lac.

6. **MFI stays OFF** (`useOffres=false`) — dropped as a correction.

## Ground truth (mM) — for per-compartment validation
| compartment | Act | Cr | Cho | Lac |
|---|---|---|---|---|
| Small Inner (70 mL) | 40.42 | 32.60 | 10.08 | 4.09 |
| Large Inner (650 mL) | 6.09 | 16.03 | 14.48 | 8.01 |
| Outer (3800 mL) | 20.01 | 16.00 | 4.81 | 0.00 |

Validate by taking **ROI medians inside each compartment** (segment_inner_phantom)
and comparing to the matching row. Outer Lac = 0 is a free negative control.

## Still open / to verify
- **Center dropout:** read the `|DC|/|k0 signal|` printout to confirm DC vs coil/B1.
- **Choline = 9 protons** (N(CH₃)₃) vs Act/Cr/Lac = 3. Verify `make_rosette_phantom_basis`
  scaleFactors, else Cho reads ~3× off.
- **TE/basis:** if the warning fires, rebuild the basis at the data TE.
