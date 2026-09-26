# EEMPrep

Modular EEM (Excitation-Emission Matrix) preprocessing library for biomedical spectroscopy.

## Installation

```bash
pip install eemprep
```

For development:

```bash
git clone https://github.com/alexgaarciia/EEMPrep.git
cd EEMPrep
pip install -e ".[dev]"
```

## Usage

```python
from eemprep.preprocessing import subtract_blank, normalize_eem, remove_rayleigh

# eem: (n_ex, n_em) matrix; ex_axis / em_axis: wavelengths in nm
eem = subtract_blank(eem, blank_scans)                    # remove solvent background
eem = normalize_eem(eem, raman_scans, raman_em_axis)      # convert to Raman Units (R.U.)
eem = remove_rayleigh(eem, ex_axis, em_axis, half_width=25)  # mask Rayleigh scatter with NaN
```

### Available steps

| Function | Module | Description |
|---|---|---|
| `subtract_blank` | `blank_subtraction` | Subtracts a blank (or the mean of several blanks) to remove background fluorescence. |
| `normalize_eem` | `raman_normalization` | Divides by the water Raman peak area (Ex = 350 nm, Em 371–428 nm). |
| `remove_rayleigh` | `scatter_removal` | Masks the region em < ex and the 1st/2nd order Rayleigh scatter bands with NaN. |

## Requirements

- Python >= 3.9
- NumPy >= 2.0

## License

MIT
