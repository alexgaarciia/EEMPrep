import numpy as np


def subtract_blank(eem, blank, clip_negative=True):
    """Subtract a blank (solvent) measurement from an EEM to remove background fluorescence.

    Parameters
    ----------
    eem : np.ndarray
        EEM matrix of shape (n_emission, n_excitation).
    blank : np.ndarray
        Blank measurement matrix of the same shape as `eem`.

    Returns
    -------
    np.ndarray
        Background-corrected EEM of the same shape as `eem`.

    Raises
    ------
    ValueError
        If `eem` and `blank` do not have the same shape.
    """
    if eem.shape != blank.shape:
        raise ValueError(f"Shape mismatch: eem {eem.shape} vs blank {blank.shape}")

    result = eem - blank
    if clip_negative:
        result = np.maximum(result, 0)

    return result
