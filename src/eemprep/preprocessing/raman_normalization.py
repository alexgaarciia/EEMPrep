import numpy as np


def normalize_eem(eem, raman_scans, raman_em_axis):
    """Normalize an EEM by the Raman scatter peak area of pure water (Ex=350 nm).

    Divides the EEM by the integrated area of the water Raman peak, making
    fluorescence intensities comparable across instruments and measurement sessions

    Parameters
    ----------
    eem : np.ndarray
        EEM matrix of shape (n_ex, n_em).
    raman_scans : list of np.ndarray
        List of water Raman scans, each of shape (n_replicates, n_em).
        Multiple scans are averaged to reduce noise.
    raman_em_axis : np.ndarray
        Emission wavelength axis for the Raman scans, shape (n_em,).

    Returns
    -------
    np.ndarray
        Raman-normalized EEM in Raman Units (R.U.), same shape as `eem`.

    Raises
    ------
    ValueError
        If the computed Raman peak area is zero or negative.
    """
    # Average replicates within each scan, then average across scans
    raman_scan_means = [scan.mean(axis=0) for scan in raman_scans]
    raman_mean = np.mean(raman_scan_means, axis=0)

    # Water Raman peak at Ex=350 nm falls between 371–428 nm (standard integration window)
    mask = (raman_em_axis >= 371) & (raman_em_axis <= 428)
    raman_roi = raman_mean[mask]
    raman_area = np.trapezoid(raman_roi, raman_em_axis[mask])

    if raman_area <= 0:
        raise ValueError(
            f"Raman peak area is {raman_area:.4f}; check that raman_scans cover 371–428 nm"
        )

    return eem / raman_area
    