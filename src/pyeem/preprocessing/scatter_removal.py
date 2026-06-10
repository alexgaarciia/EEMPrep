import numpy as np


def remove_rayleigh(eem, ex_axis, em_axis, half_width=25):
    """Zero out Rayleigh scatter bands (1st and 2nd order) from an EEM.

    Sets to zero the physically impossible region (em < ex) and a symmetric
    band around each scatter diagonal.

    Parameters
    ----------
    eem : np.ndarray
        EEM array of shape (n_samples, n_ex, n_em) or (n_ex, n_em).
    ex_axis : np.ndarray
        Excitation wavelengths, shape (n_ex,).
    em_axis : np.ndarray
        Emission wavelengths, shape (n_em,).
    half_width : int or float, optional
        Half-width in nm of the band masked around each scatter diagonal.
        Default is 25.

    Returns
    -------
    np.ndarray
        Copy of `eem` with scatter bands zeroed out.

    Raises
    ------
    ValueError
        If axis lengths do not match the corresponding EEM dimensions.
    """
    eem = eem.copy().astype(float)

    if eem.ndim == 2:
        if eem.shape != (len(ex_axis), len(em_axis)):
            raise ValueError(
                f"2D eem shape {eem.shape} does not match (n_ex={len(ex_axis)}, n_em={len(em_axis)})"
            )
    elif eem.ndim == 3:
        if eem.shape[1:] != (len(ex_axis), len(em_axis)):
            raise ValueError(
                f"3D eem shape {eem.shape} does not match (*, n_ex={len(ex_axis)}, n_em={len(em_axis)})"
            )
    else:
        raise ValueError(f"eem must be 2D or 3D, got {eem.ndim}D")

    # Grid of wavelength differences
    delta = em_axis[np.newaxis, :] - ex_axis[:, np.newaxis]

    # Emission below excitation is physically impossible (violates energy conservation)
    below_ex = delta < 0

    # 2nd-order Rayleigh: grating diffracts excitation light at exactly 2x the wavelength
    second_order = np.abs(em_axis[np.newaxis, :] - 2 * ex_axis[:, np.newaxis]) <= half_width

    mask = below_ex | (np.abs(delta) <= half_width) | second_order

    eem[..., mask] = np.nan

    return eem
