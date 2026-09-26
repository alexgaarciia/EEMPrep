import numpy as np


def subtract_blank(eem, blank_scans, clip_negative=True):
    """Subtract a blank (solvent) measurement from an EEM to remove background fluorescence.

    Parameters
    ----------
    eem : np.ndarray
        EEM matrix of shape (n_emission, n_excitation).
    blank_scans : np.ndarray
        Single blank of shape (n_emission, n_excitation) or a stack of blanks
        of shape (n_scans, n_emission, n_excitation). Multiple blanks are averaged.

    Returns
    -------
    np.ndarray
        Background-corrected EEM of the same shape as `eem`.

    Raises
    ------
    ValueError
        If `eem` and the blank do not have the same shape.
    """
    blank_scans = np.asarray(blank_scans)

    if blank_scans.ndim == 3:
        blank = np.mean(blank_scans, axis=0)
    else:
        blank = blank_scans

    if eem.shape[-2:] != blank.shape:
        raise ValueError(f"Shape mismatch: eem {eem.shape} vs blank {blank.shape}")

    result = eem - blank
    if clip_negative:
        result = np.maximum(result, 0)

    return result
