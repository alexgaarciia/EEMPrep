# EEMPrep

**EEMPrep** (*Excitation-Emission Matrix Preprocessing*) is a modular Python library for preprocessing **fluorescence EEMs** (Excitation-Emission Matrices), with a focus on biomedical spectroscopy.

It provides simple, NumPy-based building blocks to turn raw EEMs into clean, comparable data: **blank subtraction**, **Raman normalization** and **Rayleigh scatter removal**. Each step is a single function that can be used on its own or chained with the others.

<p align="center">
  <img src="https://raw.githubusercontent.com/alexgaarciia/EEMPrep/main/docs/images/overview.png" alt="EEMPrep overview" width="900">
</p>

> 🚧 **Status: early development.** The API may still change between versions.

---

## Table of contents

- [Why EEMPrep?](#why-eemprep)
- [Installation](#installation)
- [Input data format](#input-data-format)
- [Quick start](#quick-start)
- [Preprocessing steps explained](#preprocessing-steps-explained)
  - [1. Blank subtraction](#1-blank-subtraction)
  - [2. Raman normalization](#2-raman-normalization)
  - [3. Rayleigh scatter removal](#3-rayleigh-scatter-removal)
- [Recommended order](#recommended-order)
- [Repository structure](#repository-structure)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

---

## Why EEMPrep?

An **Excitation-Emission Matrix** is a 2D fluorescence map: the sample is excited at a series of wavelengths and, for each one, the full emission spectrum is recorded. The result is a "fingerprint" of the fluorophores in the sample (proteins, NADH, flavins, porphyrins, humic-like substances…), widely used in biomedical, environmental and food analysis.

Raw EEMs, however, contain signals that do not come from the sample's fluorescence:

- **Background fluorescence** from the solvent, cuvette or buffer.
- **Instrument-dependent intensities**: the same sample gives different values on different instruments, days or lamp ages.
- **Scatter bands**: Rayleigh scatter (emission = excitation) and its 2nd-order harmonic (emission = 2 × excitation) appear as intense diagonal lines that dominate the matrix and break multivariate models such as PARAFAC.

EEMPrep collects the standard corrections for these problems in one place with a consistent interface:

- Every function works on **NumPy arrays**, either a single EEM `(n_ex, n_em)` or, where indicated, a stack of EEMs `(n_samples, n_ex, n_em)`.
- Every function **returns a new array** and never modifies the input.
- Steps are **independent**, so you can use only the ones you need.

---

## Installation

EEMPrep requires **Python ≥ 3.9** and depends only on `numpy` (≥ 2.0).

### Option A: Install from PyPI

```bash
pip install eemprep
```

To install the latest development version directly from GitHub instead:

```bash
pip install git+https://github.com/alexgaarciia/EEMPrep.git
```

### Option B: Clone and install in development (editable) mode

This is the recommended option if you want to modify the code. Changes to the source files take effect immediately without reinstalling.

```bash
git clone https://github.com/alexgaarciia/EEMPrep.git
cd EEMPrep

# (optional but recommended) create a virtual environment
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

pip install -e .
```

To also install the development tools (`pytest`, `black`, `ruff`):

```bash
pip install -e ".[dev]"
```

### Check the installation

```python
import eemprep
print(eemprep.__version__)
```

---

## Input data format

All functions follow the same conventions:

| Object | Type | Shape | Description |
|---|---|---|---|
| `eem` | `numpy.ndarray` | `(n_ex, n_em)` | A single EEM. **Rows are excitation wavelengths, columns are emission wavelengths.** |
| `eem` (stack) | `numpy.ndarray` | `(n_samples, n_ex, n_em)` | Several EEMs measured on the same wavelength grid. |
| `ex_axis` | `numpy.ndarray` | `(n_ex,)` | Excitation wavelengths, in **nm**. |
| `em_axis` | `numpy.ndarray` | `(n_em,)` | Emission wavelengths, in **nm**. |

Example of loading an EEM from a CSV where the first column holds the excitation wavelengths and the header holds the emission wavelengths:

```python
import numpy as np
import pandas as pd

df = pd.read_csv("sample_eem.csv", index_col=0)

eem = df.to_numpy(dtype=float)                 # (n_ex, n_em)
ex_axis = df.index.astype(float).to_numpy()    # (n_ex,)
em_axis = df.columns.astype(float).to_numpy()  # (n_em,)
```

> **Note:** some instruments export EEMs with **emission in rows and excitation in columns**. In that case transpose the matrix (`eem = eem.T`) and swap the axes so that rows correspond to excitation.

---

## Quick start

```python
from eemprep.preprocessing import subtract_blank, normalize_eem, remove_rayleigh

# 1. Remove the solvent background (several blank scans are averaged)
eem = subtract_blank(eem, blank_scans)

# 2. Convert intensities to Raman Units (R.U.) using a water Raman scan
eem = normalize_eem(eem, raman_scans, raman_em_axis)

# 3. Mask the Rayleigh scatter bands with NaN
eem = remove_rayleigh(eem, ex_axis, em_axis, half_width=25)
```

---

## Preprocessing steps explained

This section explains **what** each step does, **why** it is useful and **when** to use it.

### 1. Blank subtraction

**Module:** `eemprep.preprocessing.blank_subtraction` · **Function:** `subtract_blank(eem, blank_scans, clip_negative=True)`

A **blank** is an EEM of the solvent or buffer alone (for example ultrapure water or PBS), measured under the same conditions as the samples. Subtracting it removes the background fluorescence of the medium and cuvette, as well as part of the Raman scatter of the solvent.

- `blank_scans`: a single blank `(n_ex, n_em)` or a stack of blanks `(n_scans, n_ex, n_em)`. When several blanks are given, they are **averaged** before subtraction, which reduces noise.
- `clip_negative`: if `True` (default), negative values produced by the subtraction are set to 0, since fluorescence intensity cannot be negative.

`eem` can be a single EEM or a stack of EEMs; in both cases the same blank is subtracted from every matrix. The blank must have the same wavelength grid as the samples, otherwise a `ValueError` is raised.

```python
from eemprep.preprocessing import subtract_blank

eem_corr = subtract_blank(eem, blank_scans)                          # blank_scans: (3, n_ex, n_em)
eem_raw_diff = subtract_blank(eem, blank_scans, clip_negative=False)  # keep negative values
```

### 2. Raman normalization

**Module:** `eemprep.preprocessing.raman_normalization` · **Function:** `normalize_eem(eem, raman_scans, raman_em_axis)`

Fluorescence intensities are given in arbitrary units that depend on the instrument, the lamp and the detector settings. To make EEMs comparable across **instruments and measurement sessions**, they are divided by the area of the **water Raman peak**, a stable reference signal. The result is expressed in **Raman Units (R.U.)**.

How it works:

1. Each water Raman scan (emission spectrum at **Ex = 350 nm**) is averaged over its replicates, and then all scans are averaged together.
2. The averaged spectrum is integrated (trapezoidal rule) over the Raman peak window, **Em = 371–428 nm**.
3. The EEM is divided by that area.

Parameters:

- `raman_scans`: a **list** of arrays, each of shape `(n_replicates, n_em)`. For example, one array per Raman measurement session.
- `raman_em_axis`: emission wavelengths of the Raman scans, in nm. It must cover the 371–428 nm window, otherwise the computed area is zero and a `ValueError` is raised.

```python
from eemprep.preprocessing import normalize_eem

# two Raman measurements of ultrapure water, 3 replicates each
raman_scans = [raman_morning, raman_afternoon]   # each: (3, n_em_raman)
eem_ru = normalize_eem(eem, raman_scans, raman_em_axis)
```

> **Tip:** measure the water Raman scan on the **same day** and with the **same slit widths and integration time** as your samples. The normalization only corrects for differences in instrument response if the reference is acquired under identical settings.

### 3. Rayleigh scatter removal

**Module:** `eemprep.preprocessing.scatter_removal` · **Function:** `remove_rayleigh(eem, ex_axis, em_axis, half_width=25)`

Part of the excitation light is scattered elastically by the sample and reaches the detector at the **same wavelength** (1st-order Rayleigh, Em = Ex). The monochromator grating also lets it through at **twice the wavelength** (2nd-order Rayleigh, Em = 2 × Ex). Both appear as intense diagonal bands that carry no chemical information and can be orders of magnitude stronger than the fluorescence.

`remove_rayleigh` sets to **NaN**:

| Region | Condition | Reason |
|---|---|---|
| Below the 1st-order line | Em < Ex | Physically impossible: emitted photons cannot have more energy than the excitation photons. |
| 1st-order Rayleigh band | \|Em − Ex\| ≤ `half_width` | Elastic scatter at the excitation wavelength. |
| 2nd-order Rayleigh band | \|Em − 2·Ex\| ≤ `half_width` | Grating harmonic of the scattered light. |

- `half_width`: half-width in nm of the band removed around each scatter line (default 25 nm). Increase it if scatter residues remain visible at the edges of the bands; decrease it if too much fluorescence is being removed.
- Works on a single EEM `(n_ex, n_em)` or a stack `(n_samples, n_ex, n_em)`.
- The output is always a **float** copy of the input.

```python
from eemprep.preprocessing import remove_rayleigh

eem_clean = remove_rayleigh(eem, ex_axis, em_axis, half_width=20)
```

> **Why NaN and not 0?** Setting the scatter regions to 0 would tell a model that there is *no fluorescence* there, which is false: the value is simply unknown. NaN marks them as **missing values**, which is how PARAFAC and other multiway tools usually expect scatter to be handled. If a later step needs finite values, fill or interpolate the NaNs explicitly (e.g. `np.nan_to_num(eem_clean)`).

---

## Recommended order

A typical workflow for a batch of samples is:

1. **Blank subtraction**: remove the solvent background.
2. **Raman normalization**: convert to Raman Units so that samples from different sessions are comparable.
3. **Rayleigh scatter removal**: mask the scatter bands as the last step, so that NaNs do not propagate through the previous operations.

```python
from eemprep.preprocessing import subtract_blank, normalize_eem, remove_rayleigh

eems = subtract_blank(eems, blank_scans)            # eems: (n_samples, n_ex, n_em)
eems = normalize_eem(eems, raman_scans, raman_em_axis)
eems = remove_rayleigh(eems, ex_axis, em_axis)
```

---

## Repository structure

```
EEMPrep/
├── src/
│   └── eemprep/
│       ├── __init__.py
│       └── preprocessing/
│           ├── __init__.py
│           ├── blank_subtraction.py      # solvent blank subtraction
│           ├── raman_normalization.py    # normalization to Raman Units
│           └── scatter_removal.py        # 1st/2nd order Rayleigh scatter masking
├── docs/
│   └── images/                           # figures used in this README
├── pyproject.toml
├── LICENSE
└── README.md
```

---

## Roadmap

- [ ] Inner filter effect (IFE) correction
- [ ] Raman scatter removal
- [ ] Interpolation of the removed scatter regions
- [ ] Configurable preprocessing pipeline
- [ ] Plotting utilities (contour maps of EEMs)
- [ ] Unit tests
- [x] Publication on PyPI

---

## Contributing

Contributions, bug reports and suggestions are welcome. Please open an [issue](https://github.com/alexgaarciia/EEMPrep/issues) or a pull request.

For development:

```bash
pip install -e ".[dev]"
ruff check src
black src
pytest
```

---

## Citation

If you use EEMPrep in your research, please cite it:

```bibtex
@software{garcia_navarro_eemprep,
  author  = {García Navarro, Alejandro Leonardo},
  title   = {EEMPrep: modular preprocessing for fluorescence excitation-emission matrices},
  year    = {2026},
  url     = {https://github.com/alexgaarciia/EEMPrep}
}
```

---

## License

This project is licensed under the **MIT License**. See [LICENSE](https://github.com/alexgaarciia/EEMPrep/blob/main/LICENSE) for details.
