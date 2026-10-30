from typing import List, Optional, Union
import numpy as np
import pyshtools
import xarray as xa


# ==============================================================================
# Helpful functions
# ==============================================================================
def ensure_dh_grid(da: xa.DataArray) -> xa.DataArray:
    """fx: ensure_dh_grid()
    Ensures the grid satisfies Driscoll-Healy dimensions for pyshtools
    (sampling=2).
    sampling=2 requires nlon = 2 * nlat; if lat shape is nlon/2+1,
    usually it runs from 90 to -90, so sort ascending and slice off 1 lat
    (-90) to match DH2 requirements.
    """
    da = da.sortby("latitude", ascending=False)
    nlat, nlon = da.shape[-2], da.shape[-1]

    if nlon == 2 * (nlat - 1):
        da = da.isel(latitude=slice(0, nlat - 1))

    return da


def splice_and_interpolate(
    lam_da: xa.DataArray,
    glm_da: xa.DataArray,
) -> xa.DataArray:
    """fx: splice_and_interpolate()
    Splices regional model output data (lam_da) into its driving global
    model data (glm_da).
    Returns the global model data with the regional model data inserted.
    """

    # Ensure coordinates are in the right order
    lam_da = lam_da.sortby("latitude", ascending=True)
    glm_da = glm_da.sortby("latitude", ascending=True)

    # Splice local regional domain into global data
    spliced_da = glm_da.copy()
    spliced_da.loc[{"longitude": lam_da.longitude, "latitude": lam_da.latitude}] = (
        lam_da
    )

    # Interpolate any missing values
    spliced_da = spliced_da.interpolate_na(
        dim="latitude", method="linear"
    ).interpolate_na(dim="longitude", method="linear")

    return spliced_da


def mass_flux(ds):
    g = 9.80665
    cos_lat = xa.ufuncs.cos(xa.ufuncs.radians(ds.latitude))
    meridional_overturning_vertical_mass_flux = (
        (-ds.meridional_overturning_vertical_velocity * cos_lat / g)
        .rename("meridional_overturning_vertical_mass_flux")
        .assign_attrs(units="kg m-2 s-1")
    )
    zonal_overturning_vertical_mass_flux = (
        (-ds.zonal_overturning_vertical_velocity * cos_lat / g)
        .rename("zonal_overturning_vertical_mass_flux")
        .assign_attrs(units="kg m-2 s-1")
    )
    return ds.assign(
        meridional_overturning_vertical_mass_flux=meridional_overturning_vertical_mass_flux,
        zonal_overturning_vertical_mass_flux=zonal_overturning_vertical_mass_flux,
    )


# ==============================================================================
# Spectral Overturning Calculation
# ==============================================================================
def compute_overturning(
    omega: xa.DataArray,
    truncs: Optional[Union[int, List[Optional[int]]]] = None,
) -> xa.Dataset:
    """fx: compute_overturning()
    Computes the vector streamfunction "mu" and returns overturning
    components across all truncations for a 2D slice.

    truncs may be None (use the grid's natural lmax), a single int, or a
    list of ints. A list gives a "truncation" dimension; None or a single
    int gives no truncation dimension, only a scalar "truncation" coord.

    Outputs are on the extended DH grid (both poles, 0 to 360-dlon), in
    the form of an xarray.Dataset with variables for each component of
    overturning streamfunction and overturning vertical velocity.

    """

    # Enforce a Driscoll-Healy grid with sampling type 2
    # in terms of the pyshtool settings (i.e. with nlon=2*nlat)
    omega = ensure_dh_grid(omega)

    # Use PySHTOOLS to extract spherical harmonic coefficients for this dataset.
    # Specifically this returns a grid of shape (2, lmax+1, lmax+1) where lmax is
    # the maximum possible expansion degree for this dataset's resolution. The first
    # axis specifies coefficients for cosines and sines (length 2), the latter two
    # specifying the values for each "l" and "m".
    omega_coeffs = pyshtools.expand.SHExpandDH(omega.values, sampling=2)

    # Identify maximum possible expansion degree for this grid - the maximum
    # "l" in the spherical harmonics that is available at this resolution
    max_available_l = omega_coeffs.shape[1] - 1

    # Normalise truncs to a list of ints, with None meaning the natural lmax
    scalar_trunc = truncs is None or isinstance(truncs, (int, np.integer))
    if scalar_trunc:
        truncs = [truncs]
    truncs = [max_available_l if t is None else int(t) for t in truncs]

    # Define the fields we will calculate, with their attributes. The Poisson
    # solve is on the unit sphere (a = 1), so psi is psi/a and carries the
    # units of omega; omega_lam and omega_phi are physical
    omega_units = omega.attrs.get("units", "Pa s-1")
    psi_comment = (
        "Unit-sphere vector streamfunction (psi / a); multiply by the Earth "
        f"radius a in m for physical units of {omega_units} m"
    )
    var_attrs = {
        "zonal_overturning_streamfunction": {
            "units": omega_units,
            "comment": psi_comment,
        },
        "meridional_overturning_streamfunction": {
            "units": omega_units,
            "comment": psi_comment + "; positive southward (d/d colatitude)",
        },
        "zonal_overturning_vertical_velocity": {"units": omega_units},
        "meridional_overturning_vertical_velocity": {"units": omega_units},
    }
    fields = {name: [] for name in var_attrs}

    # Loop over desired truncations to do calculations at various lmax
    for trunc in truncs:
        # define lmax for the calculation which is the lmax specified in the
        # truncations field, or if that's too big, the maximum l for the grid
        trunc_lmax = min(trunc, max_available_l)

        # set up mu_coeffs grid
        mu_coeffs = omega_coeffs.copy()

        # for truncations smaller than the grid's natural lmax, set all
        # values of the coefficients beyond lmax's index to 0
        if trunc_lmax < max_available_l:
            mu_coeffs[:, trunc_lmax + 1 :, :] = 0.0

        # for the truncation's lmax, loop through each l and convert the
        # omega coefficients into mu coefficients by dividing by -l(l+1)
        for l in range(1, trunc_lmax + 1):  # skip l = 0
            mu_coeffs[:, l, : l + 1] /= -l * (l + 1)

        # re-expand mu and calculate its gradient on the spherical grid
        grad_mu = pyshtools.SHCoeffs.from_array(mu_coeffs).gradient()
        psi_lam = grad_mu.phi.data
        psi_phi = grad_mu.theta.data

        # expand coefficients to full grid for streamfunction, cutting out duplicate
        # coordinates in lon and lat
        psi_lam_coeffs = pyshtools.expand.SHExpandDH(psi_lam[:-1, :-1], sampling=2)
        psi_phi_coeffs = pyshtools.expand.SHExpandDH(psi_phi[:-1, :-1], sampling=2)

        # Calculate omega by taking divergence of psi
        # calculate gradients of components of psi
        grad_psi_lambda = pyshtools.SHCoeffs.from_array(psi_lam_coeffs).gradient()
        grad_psi_phi = pyshtools.SHCoeffs.from_array(psi_phi_coeffs).gradient()

        # Lambda component of divergence is just the same as the lambda component
        # of gradient
        omega_lam = grad_psi_lambda.phi.data

        # Phi component of divergence is more complex:
        # (1/sin th) d(sin th psi)/d th =  d psi/d th + cot(th) psi
        # th=theta=colatitude. cot(colat) = tan(lat)
        tanlat = np.tan(np.deg2rad(grad_psi_phi.theta.lats()))[:, None]
        omega_phi = grad_psi_phi.theta.data
        omega_phi = omega_phi + tanlat * psi_phi
        omega_phi[[0, -1], :] = np.nan  # undefined at poles. tan(+-pi/2)

        # Construct fields of the dataset...
        for name, field in zip(var_attrs, (psi_lam, psi_phi, omega_lam, omega_phi)):
            # Drop the repeated 360E column of the extended output grid
            fields[name].append(field[:, :-1])

    # Keep any scalar coords of the input slice (e.g. pressure, time)
    scalar_coords = {k: v for k, v in omega.coords.items() if v.ndim == 0}

    # Construct dataset
    ds = xa.Dataset(
        data_vars={
            name: (
                ("truncation", "latitude", "longitude"),
                np.stack(fields[name]),
                attrs,
            )
            for name, attrs in var_attrs.items()
        },
        coords={
            "truncation": (
                "truncation",
                truncs,
                {"units": "1", "long_name": "maximum spherical harmonic degree"},
            ),
            "latitude": grad_mu.phi.lats(),
            "longitude": grad_mu.phi.lons()[:-1],
            **scalar_coords,
        },
    )

    # Remove truncation dimension if it is unnecessary
    if scalar_trunc:
        ds = ds.squeeze("truncation")

    return ds
