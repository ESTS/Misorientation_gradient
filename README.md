# Misorientation gradient
## kamslope — misorientation gradient from EBSD orientation maps

Python implementation of the KAM-slope method used to quantify the local
misorientation gradient (Δθ/Δx, in °/µm) from square-grid EBSD scans, as a
proxy for the stored energy of plastic deformation.

## Scope

This repository contains the calculation of a single quantity: the slope of
the kernel average misorientation with respect to kernel radius, evaluated
per pixel, with cubic crystal symmetry.

It does **not** contain texture analysis, ODF calculation, orientation
binning in Euler space, grain reconstruction, or GND density inversion. For
those, use MTEX or the OIM software.

## Method

For each pixel, the kernel average misorientation KAM(n) is evaluated for
kernels of radius n = 1 … n_max in grid steps. KAM increases linearly with
kernel radius over a limited range of n; the slope of that line, expressed
in °/µm, corresponds to a component of the lattice curvature tensor and
scales with the density of geometrically necessary dislocations. Taking the
slope removes the dependence of a single KAM value on the arbitrary choice
of kernel size.

Neighbour pairs above a cut-off (5° by default) are excluded, so that
(sub)grain boundaries do not contribute to the intragranular average.

## Two caveats that affect the reported number

**1. The absolute value depends on the kernel definition.** Two kernels are
implemented: `perimeter` (neighbours at Chebyshev distance exactly n) and
`filled` (all neighbours within distance n). They do not give the same slope,
because the mean neighbour separation inside the kernel differs. For the
perimeter kernel, the mean Euclidean neighbour distance is 1.1269 · n · step,
not n · step; the returned slope is larger than the true lattice curvature by
that factor. The example scan below reproduces this exactly: an imposed
curvature of 5.000 °/µm is returned as 5.635 °/µm.

Slopes are therefore only comparable between studies that state the kernel
geometry. No implicit calibration is applied here.

**2. The intercept is not forced through the origin by default.** A free
intercept absorbs the angular noise floor of the measurement rather than
pushing it into the slope. Use `force_origin=True` to change this, and report
which was used.

Angular noise in the raw orientation data propagates into KAM(1) most
strongly. If the data have not been filtered (e.g. Kuwahara), inspect the
linearity plot before trusting the slope.

## Two routes

**A. From orientations** (`kamslope.kam`) — computes KAM from the Euler angles
of a square-grid `.ang` scan, then fits the slope.

**B. From OIM KAM exports** (`kamslope.oim`) — takes KAM maps already exported
from OIM for kernel orders n = 1, 2, 3 … and fits the slope against the
neighbour distance. This is the route used in practice, and it supports both
square and hexagonal grids.

### Grid geometry and kernel distances

The abscissa of the fit is the neighbour distance of each kernel order, and it
depends on the lattice:

| shell | square | hexagonal |
|---|---|---|
| 1st | 1 | 1 |
| 2nd | √2 ≈ 1.414 | √3 ≈ 1.732 |
| 3rd | 2 | 2 |
| 4th | √5 ≈ 2.236 | √7 ≈ 2.646 |

Only the first and third shells coincide. Shell distances and multiplicities
are enumerated from the lattice vectors, not tabulated, so they stay correct
beyond the third neighbour.

The geometry is selected with `--grid {auto,square,hexagonal}`. `auto` infers
it from the exported coordinates and is the default. An explicit choice is
still checked against the coordinates and is rejected if the two disagree;
`--force-grid` downgrades that to a warning. This asymmetry is deliberate:
declaring a hexagonal scan to be square biases the slope upward by a few per
cent, in one direction, with no symptom in any downstream quantity, so the
declaration has to fail loudly rather than be trusted.

```bash
# geometry inferred from the coordinates
python scripts/slope_from_oim_export.py KAM1.txt KAM2.txt KAM3.txt \
       --ecdf-out slope_ecdf.txt

# geometry and step declared explicitly, still cross-checked
python scripts/slope_from_oim_export.py KAM1.txt KAM2.txt KAM3.txt \
       --grid hexagonal --step 1.0 --kernel perimeter
```

`--kernel perimeter` assumes OIM's perimeter-only option was active, so that
KAM(n) samples the n-th shell alone. `--kernel cumulative` is for the case
where KAM(n) averages over all shells up to n; the abscissa then becomes the
occupancy-weighted mean distance over the kernel (1, 1.366, 1.577 for a
hexagonal grid). On the same data the two differ by close to a factor of two,
so the setting used must be recorded.

Points at the ceiling of the OIM KAM definition (5° by default) and points at
exactly 0.000° are excluded from the fit: they are censored values rather than
measurements, and including them leaves the mean almost unchanged while
inflating the standard deviation by roughly a quarter.

## Installation

```bash
git clone https://github.com/<user>/kamslope.git
cd kamslope
pip install -e .
```

Requires Python ≥ 3.10, NumPy and Matplotlib. No compiled dependencies.

## Usage

```bash
# runs on the bundled example scan and writes a figure + summary JSON
python scripts/run_kam_slope.py data/example_scan.ang --out figures/

# no argument: uses a synthetic scan of known curvature
python scripts/run_kam_slope.py
```

As a library:

```python
from kamslope import read_ang, kam_slope

scan = read_ang("my_scan.ang")
res = kam_slope(scan.euler, scan.step_um, n_max=4,
                kernel="perimeter", threshold_deg=5.0,
                mask=scan.ci_mask(0.1))
print(res.summary())          # mean, median, std in deg/um
res.slope                     # per-pixel map
res.r_squared                 # per-pixel fit quality
```

## Validation

```bash
pytest tests -q
```

The suite checks the parts that fail silently rather than loudly:

- the cubic group is generated, not assumed, and contains 24 proper rotations;
- the closed-form disorientation angle agrees with an explicit minimisation
  over all 576 symmetry pairs;
- an orientation and its 24 symmetric equivalents return zero disorientation;
- no cubic–cubic pair exceeds the 62.8° Mackenzie limit;
- on a synthetic field of prescribed curvature, KAM is linear in kernel radius
  (R² > 0.999), the slope scales with the imposed curvature, and it is
  invariant under a rigid rotation of the sample frame;
- a strain-free scan returns exactly zero.

The rotation-invariance test is the one that matters most: lattice curvature
is a property of the crystal, so any dependence of the result on the sample
reference frame is a bug in the symmetry handling.

## Data

`data/example_scan.ang` is a synthetic 120 × 120 scan at 0.1 µm step with a
uniform curvature of 5 °/µm and 0.05° angular noise. It exists so the
pipeline can be verified against a known answer. It is not experimental data.

Experimental EBSD datasets are archived separately; see `data/README.md`.

## Citation

See `CITATION.cff`. If this code contributes to published work, please cite
both the software record and the paper in which the method is applied.

## Licence

MIT. See `LICENSE`.

