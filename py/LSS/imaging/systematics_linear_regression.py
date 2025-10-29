"""Taken from https://github.com/domichbt/alaeboss/tree/d14124894cf6738ae6ec11a948e2e14b3a3e0a28/src/alaeboss by fusing `linear_regressor.py` and `produce_imweights.py`."""

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

import logging
from collections.abc import Sequence
from functools import partial
from pathlib import Path
from time import time

import fitsio
import healpy as hp
import jax
import jax.numpy as jnp
import LSS.common_tools as common
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from jax.scipy.optimize import minimize as jaxminimize
from jax.typing import ArrayLike


@jax.jit
def _my_bincount(idx, accumutalor, weights):
    return accumutalor.at[idx].add(weights)[1:-1]


class LinearRegressor:
    """Object that contains the data and randoms information, and can perform à la eBOSS linear regression as well as export weights."""

    constant: str = "constant"
    """Name of the constant contribution in the density model."""
    bin_margin: float = 1e-7
    """Extra space at the beginning and end of edges to avoid null-sized bins"""

    logger = logging.getLogger(__qualname__)
    """Logger object for the class."""
    logger.setLevel(logging.DEBUG)
    _stream_handler = logging.StreamHandler()
    _formatter = logging.Formatter(
        "%(asctime)s [%(name)s] %(levelname)s | %(message)s", "%Y-%m-%d %H:%M:%S"
    )
    _stream_handler.setFormatter(_formatter)
    logger.addHandler(_stream_handler)

    def __init__(
        self,
        data_weights: ArrayLike,
        random_weights: ArrayLike,
        templates: dict[str, tuple[ArrayLike, ArrayLike]],
        loglevel: str = "INFO",
    ):
        """
        Initialize a LinearRegressor instance for imaging systematics, by giving information on the data and randoms' properties.

        The regression can then be run by calling the ``cut_outliers``, ``prepare`` and ``regress`` methods in that order. Once the regression has succeeded, weights can then be obtained by calling the ``weight_model`` method on the ``coefficients`` attribute.

        Parameters
        ----------
        data_weights : ArrayLike
            External, non-imaging weights for each object of the data catalog (FKP, completeness...)
        random_weights : ArrayLike
            External, non-imaging weights for each object of the randoms catalog (FKP, completeness...)
        templates : dict[str, tuple[ArrayLike, ArrayLike]]
            Dictionnary containing the imaging template values for each objects of the data and randoms catalogs, under the form of {template_name: (data_values, random_values)}
        loglevel : str
            Minimal logging level, set by default to INFO. Additional output can be obtained from DEBUG.
        """
        # Build the systematics matrices
        template_values_data = jnp.vstack(
            [jnp.array(template[0]) for template in templates.values()]
        )
        template_values_randoms = jnp.vstack(
            [jnp.array(template[1]) for template in templates.values()]
        )

        self._initial_setup(
            data_weights=data_weights,
            random_weights=random_weights,
            template_values_data=template_values_data,
            template_values_randoms=template_values_randoms,
            template_names=list(templates.keys()),
            loglevel=loglevel,
        )

    def _initial_setup(
        self,
        data_weights: ArrayLike,
        random_weights: ArrayLike,
        template_values_data: ArrayLike,
        template_values_randoms: ArrayLike,
        template_names: list[str],
        loglevel: str = "INFO",
    ):
        self.logger.setLevel(loglevel)

        # Test that data types are consistent with current setting
        use_float64 = jax.config.jax_enable_x64
        if not use_float64:
            self.logger.warning(
                "``jax_enable_x64`` is set to False; this means that computations will be capped to float32. The regression will most likely fail and at least return unexpected results. Are you sure you meant to do this? Otherwise use ``jax.config.update('jax_enable_x64', True)"
            )

        self.template_names = template_names
        self.logger.info("Setting up %i templates", len(self.template_names))

        self.coefficients = jnp.full(
            (len(self.template_names) + 1), fill_value=jnp.nan, dtype=float
        )  # Best fit coefficients, empty for now

        # Templates values will be stored as numpy arrays to avoid for loops, keep a private translation between index and template name
        self._template_name_to_idx = {
            name: idx for idx, name in enumerate(self.template_names)
        }
        self._idx_to_template_name = dict(enumerate(self.template_names))

        # Find bad values and build mask with the shape of the data/the randoms
        self.good_values_data = jnp.invert(
            jnp.any(
                jnp.isinf(template_values_data) | jnp.isnan(template_values_data),
                axis=0,
            )
        )
        self.good_values_randoms = jnp.invert(
            jnp.any(
                jnp.isinf(template_values_randoms) | jnp.isnan(template_values_randoms),
                axis=0,
            )
        )

        # Mask the data and the randoms
        self.data = data_weights[self.good_values_data]
        self.randoms = random_weights[self.good_values_randoms]
        self.normalization = jnp.sum(self.randoms) / jnp.sum(self.data)

        self.template_values_data = template_values_data[:, self.good_values_data]
        self.template_values_randoms = template_values_randoms[
            :, self.good_values_randoms
        ]
        self._template_values_data_for_export = self.template_values_data.copy()

        self.logger.info("Removed bad systematics values:")
        self.logger.info(
            "Data:    %i -> %i", self.good_values_data.size, self.data.size
        )
        self.logger.info(
            "Randoms: %i -> %i", self.good_values_randoms.size, self.randoms.size
        )
        self.logger.debug("Current normalization factor is %f", self.normalization)

    @classmethod
    def from_stacked_templates(
        cls,
        data_weights: ArrayLike,
        random_weights: ArrayLike,
        template_values_data: ArrayLike,
        template_values_randoms: ArrayLike,
        template_names: list[str],
        loglevel: str = "INFO",
    ) -> Self:
        """
        Initialize a LinearRegressor instance for imaging systematics, by giving information on the data and randoms' properties.

        The regression can then be run by calling the ``cut_outliers``, ``prepare`` and ``regress`` methods in that order. Once the regression has succeeded, weights can then be obtained by calling the ``weight_model`` method on the ``coefficients`` attribute.

        Parameters
        ----------
        data_weights : ArrayLike
            External, non-imaging weights for each object of the data catalog (FKP, completeness...)
        random_weights : ArrayLike
            External, non-imaging weights for each object of the randoms catalog (FKP, completeness...)
        template_values_data : ArrayLike
            Array containing the template values for the data, for each template (shape = (# of templates, length of data))
        template_values_randoms : ArrayLike
            Array containing the template values for the randoms, for each template (shape = (# of templates, length of randoms))
        template_names : list[str]
            List of the names of the templates, in the same order as in ``template_values_data`` and  ``template_values_randoms``
        loglevel : str
            Minimal logging level, set by default to INFO. Additional output can be obtained from DEBUG

        Returns
        -------
        LinearRegressor
            LinearRegressor instance, as would have been initialized with templates written in dictionaries (default ``__init__``).
        """
        reg = cls.__new__(
            cls,
            data_weights=data_weights,
            random_weights=random_weights,
            template_values_data=template_values_data,
            template_values_randoms=template_values_randoms,
            template_names=template_names,
            loglevel=loglevel,
        )
        reg._initial_setup(
            data_weights=data_weights,
            random_weights=random_weights,
            template_values_data=template_values_data,
            template_values_randoms=template_values_randoms,
            template_names=template_names,
            loglevel=loglevel,
        )
        return reg

    def save(self, filepath: str):
        """
        Save the current instance data to disk. No regression or bin related information is saved. After loading, do not repeat ``cut_outliers`` as they are already removed.

        filepath : str
            Path to the file where data should be saved.
        """
        np.savez(
            file=filepath,
            data_weights=self.data,
            random_weights=self.randoms,
            template_values_data=self.template_values_data,
            template_values_randoms=self.template_values_randoms,
            template_names=np.array(self.template_names),
        )
        self.logger.info("Saved to %s", filepath)

    @classmethod
    def load(cls, filepath: str, loglevel: str = "INFO") -> Self:
        """
        Create a Linear regressor instance from disk.

        Parameters
        ----------
        filepath : str
            Path from which to load the file.
        loglevel : str, optional
            Logging level for the new ``LinearRegressor`` instance, by default "INFO".

        Returns
        -------
        Self
            ``LinearRegressor`` instance loaded from disk.
        """
        from_disk = jnp.load(filepath)
        return cls.from_stacked_templates(
            data_weights=from_disk["data_weights"],
            random_weights=from_disk["random_weights"],
            template_values_data=from_disk["template_values_data"],
            template_values_randoms=from_disk["template_values_randoms"],
            template_names=list(from_disk["template_names"]),
            loglevel=loglevel,
        )

    def cut_outliers(self, tail: float = 0.5):
        r"""
        For each template, remove the data and random that possess values in the extremal tails of the distribution. In total, ``tail``\% of the values are removed for each template.

        Parameters
        ----------
        tail : float
            Percentage of the extremal values to remove for each template.
        """
        lower_tail = jnp.percentile(self.template_values_data, tail / 2, axis=1)
        upper_tail = jnp.percentile(self.template_values_data, 100 - tail / 2, axis=1)

        self.extremal_data_mask = jnp.invert(
            jnp.any(
                (self.template_values_data.T < lower_tail).T
                | (self.template_values_data.T > upper_tail).T,
                axis=0,
            )
        )
        extremal_random_mask = jnp.invert(
            jnp.any(
                (self.template_values_randoms.T < lower_tail).T
                | (self.template_values_randoms.T > upper_tail).T,
                axis=0,
            )
        )

        # Apply cuts
        self.data = self.data[self.extremal_data_mask]
        self.randoms = self.randoms[extremal_random_mask]
        self.normalization = jnp.sum(self.randoms) / jnp.sum(self.data)

        self.template_values_data = self.template_values_data[
            :, self.extremal_data_mask
        ]
        self.template_values_randoms = self.template_values_randoms[
            :, extremal_random_mask
        ]

        self.logger.info("Removed extremal systematics values (%f %%):", tail)
        self.logger.info(
            "Data:    %i -> %i", self.extremal_data_mask.size, self.data.size
        )
        self.logger.info(
            "Randoms: %i -> %i", extremal_random_mask.size, self.randoms.size
        )
        self.logger.debug("Current normalization factor is %f", self.normalization)

    def prepare(self, nbins: int):
        """
        Define bins for the systematics values and setup the corresponding random's binned counts.

        Parameters
        ----------
        nbins : int
            Number of bins for the imaging systematics values.
        """
        self.logger.info("Setting up %i bins", nbins)
        self.nbins = nbins
        minimal_dtype = np.min_scalar_type(
            int(nbins) + 1
        )  # save some space (most likely int8 instead of int64)
        self.edges = jnp.linspace(
            jnp.min(self.template_values_data, axis=1) - self.bin_margin,
            jnp.max(self.template_values_data, axis=1) + self.bin_margin,
            self.nbins + 1,
            axis=1,
        )

        # useful: renormalize data and randoms to [0, 1] based on edges
        self.template_values_data_normalized = (
            (self.template_values_data.T - self.edges[:, 0])
            / (self.edges[:, -1] - self.edges[:, 0])
        ).T
        template_values_randoms_normalized = (
            (self.template_values_randoms.T - self.edges[:, 0])
            / (self.edges[:, -1] - self.edges[:, 0])
        ).T
        self.logger.debug("Renormalized data and randoms")

        # random bin counts are fixed here, as they do not depend on regression coefficients
        # Equivalent to digitize but faster  # shape (N_sys, len(randoms))
        self.randoms_digitized = (
            jnp.floor(template_values_randoms_normalized * nbins).astype(minimal_dtype)
            + 1
            - (self.template_values_randoms.T == self.edges[:, -1]).T
        )
        self.logger.debug("Digitized randoms")

        # no axis-dependent implementation of bincount or histogram so have to initialize and fill array
        # jax.lax.map because the randoms array can get too big for JAX jax.vmap implementation otherwise
        self.randoms_binned = jax.lax.map(
            f=partial(
                _my_bincount,
                accumutalor=jnp.zeros(self.nbins + 2, dtype=float),
                weights=self.randoms,
            ),
            xs=self.randoms_digitized,
        )
        self.logger.debug("Binned randoms")

        # also digitize the data for later use (can simply bincount with updated model weights)
        self.data_digitized = (
            jnp.floor(self.template_values_data_normalized * nbins).astype(
                minimal_dtype
            )
            + 1
            - (self.template_values_data.T == self.edges[:, -1]).T
        )  # Equivalent to digitize but faster
        self.logger.debug("Digitized data")
        self.data_binned_noweights = jax.vmap(
            partial(
                _my_bincount,
                accumutalor=jnp.zeros(self.nbins + 2, dtype=float),
                weights=self.data,
            )
        )(self.data_digitized)
        self.logger.debug("Binned data")

        # and fix the errors on each bin
        error_everywhere = self.normalization * jnp.sqrt(
            self.data_binned_noweights / self.randoms_binned**2
            + self.data_binned_noweights**2 / self.randoms_binned**3
        )
        self.error = jnp.where(
            (self.data_binned_noweights == 0), 1e10, error_everywhere
        )
        self.logger.debug("Computed error on histogram ratio")
        self.initial_chi2 = jnp.sum(
            (self.normalization * self.data_binned_noweights / self.randoms_binned - 1)
            ** 2
            / self.error**2
        )
        self.logger.debug("Computed initial, correction-less chi2")

    def weight_model(self, coefficients: ArrayLike):
        """
        Return the weights computed from the model for the data. The returned array has the same shape as the good values of the data.

        Parameters
        ----------
        coefficients : ArrayLike
            Parameters of the regression to weight each template, as well as the constant offset at index 0.
        """
        return 1 / (
            1
            + coefficients[0]
            + jnp.dot(coefficients[1:], self.template_values_data_normalized)
        )

    def chi2(self, coefficients: ArrayLike) -> float:
        """
        Return the total cost function, defined for each template as the square of the ratio (weighted data histogram) to (weighted random histogram) minus one.

        Coefficients are used for the weight computation for the data.
        Use Poisson error of the ratio of weighted data (using coefficients) and weighted randoms.

        Parameters
        ----------
        coefficients : ArrayLike
            Parameters of the regression to weight each template, as well as the constant offset at index 0.

        Returns
        -------
        float
            chisquare value for the given coefficients.
        """
        data_binned = jax.vmap(
            partial(
                _my_bincount,
                accumutalor=jnp.zeros(self.nbins + 2, dtype=float),
                weights=self.data * self.weight_model(coefficients),
            )
        )(self.data_digitized)
        # error_everywhere = self.normalization * \
        #     jnp.sqrt(data_binned / self.randoms_binned**2 + data_binned**2 / self.randoms_binned**3) # model dependent error
        # Compute the chisquare over actual imaging templates (not the constant, which is mostly useful for the weight model)
        return jnp.sum(
            (self.normalization * data_binned / self.randoms_binned - 1) ** 2
            / self.error**2
        )

    def regress(self, guess: ArrayLike | None = None) -> dict[str, float] | None:
        """
        Find optimal coefficients by minimizing the cost function. Can provide an initial guess ``guess``, otherwise all ones will be used.

        Parameters
        ----------
        guess : ArrayLike | None
            Initial guess for the regression coefficients, of shape (1 + number of templates). Default is ``None``.

        Returns
        -------
        dict[str, float] | None
            Dictionnary of the fit coefficient for each template. The offset is referred to as 'constant'. If the minimization fails, ``None`` is returned.
        """
        if guess is None:
            guess = jnp.ones(shape=(len(self.template_names) + 1), dtype=float)

        self.logger.info("Starting minimization with initial guess %s", str(guess))
        self.logger.info(
            "Without photometry model, the chisquare is %f", self.initial_chi2
        )

        res = jaxminimize(self.chi2, guess, method="BFGS")
        if not res.success:
            self.logger.error("The optimization failed with status %i", res.status)
            return None
        else:
            self.logger.info("Minimization succeeded, chisquare = %f", self.chi2(res.x))
            self.coefficients = res.x
            return dict(zip([self.constant, *self.template_names], res.x, strict=True))

    def regress_minuit(self, guess: ArrayLike | None = None) -> dict[str, float] | None:
        """
        Find optimal coefficients by minimizing the cost function using iminuit. Can provide an initial guess ``guess``, otherwise all ones will be used.

        Parameters
        ----------
        guess : ArrayLike | None
            Initial guess for the regression coefficients, of shape (1 + number of templates). Default is ``None``.

        Returns
        -------
        dict[str, float] | None
            Dictionnary of the fit coefficient for each template. The offset is referred to as 'constant'. If the minimization fails, ``None`` is returned.
        """
        from iminuit import Minuit

        if guess is None:
            guess = jnp.zeros(shape=(len(self.template_names) + 1), dtype=float)

        self.logger.info("Starting minimization with initial guess %s", str(guess))
        self.logger.info(
            "Without photometry model, the chisquare is %f", self.initial_chi2
        )

        m = Minuit(jax.jit(self.chi2), guess, grad=jax.jit(jax.grad(self.chi2)))
        m.errordef = Minuit.LEAST_SQUARES
        m.limits = [(None, None)] * len(guess)
        m.errors = [0.1] * len(guess)

        m.migrad()
        if not m.valid:
            self.logger.error("The optimization failed with status\n%s", m.fmin)
            return None
        else:
            self.logger.info(
                "Minimization succeeded, chisquare = %f\n%s", m.fval, m.fmin
            )
            self.coefficients = jnp.array(m.values)
            return dict(
                zip([self.constant, *self.template_names], list(m.values), strict=True)
            )

    def export_weights(self) -> ArrayLike:
        """
        Return the imaging systematics weights using the parameters computed from the regression.

        Notes
        -----
        The weights are also computed for data with extremal values of the systematics that were ignored in the regression.
        """
        assert (~jnp.isnan(self.coefficients)).any()
        _template_values_data_for_export_normalized = (
            (self._template_values_data_for_export.T - self.edges[:, 0])
            / (self.edges[:, -1] - self.edges[:, 0])
        ).T
        return 1 / (
            1
            + self.coefficients[0]
            + jnp.dot(
                self.coefficients[1:], _template_values_data_for_export_normalized
            )
        )

    def mask(self) -> ArrayLike:
        """
        Return a boolean mask linking the original data catalog to the internal one (good and non-extremal values only).

        Returns
        -------
        ArrayLike
            Boolean mask with 0 for bad data or extremal systematics values and 1 for the rest.
        """
        if hasattr(self, "extremal_data_mask"):
            return self.good_values_data.at[self.good_values_data].set(
                self.extremal_data_mask
            )
        else:
            return self.good_values_data

    # Additional methods to match original ``Syst`` class
    def get_subsample(self, subdata_mask: ArrayLike) -> Self:
        """
        Create a new LinearRegressor instance from a subsample of the data.

        Parameters
        ----------
        subdata_mask : ArrayLike
            Mask to apply on the data and templates. Note that these should have had bad values and outliers cut out, so the shape may be different from initial data fed to the initializer. Correspondance between original and current shape can be obtained from ``LinearRegressor.mask()``.
        """
        data_weights = self.data[subdata_mask]
        random_weights = self.randoms
        templates = {
            name: (
                self.template_values_data[self._template_name_to_idx[name]][
                    subdata_mask
                ],
                self.template_values_randoms[self._template_name_to_idx[name]],
            )
            for name in self.template_names
        }
        new_regressor = self.__class__(
            data_weights=data_weights,
            random_weights=random_weights,
            templates=templates,
            loglevel=self.logger.level,
        )
        # do not cut outliers again
        new_regressor.prepare(nbins=self.nbins)
        return new_regressor

    def plot_overdensity(
        self,
        coefficients: ArrayLike | None = None,
        ylim: tuple[float, float] = [0.75, 1.25],
        nbinsh: int = 50,
        title: str | None = None,
    ):
        """
        Create a subplot for each template, and plot the normalized histogram values before and after the regression. The initial distribution of the template on the data is overlaid on the bottom of the plot.

        Parameters
        ----------
        coefficients : ArrayLike | None
            parameters for the weight model; if None, will default to the regression result.
        ylim : tuple[float, float]
            Limits for the y-axis
        nbinsh : int
            Number of bins for the overlaid histograms
        title : str
            Optional title
        """
        fig, axes = plt.subplots(
            1,
            len(self.template_names),
            sharey=True,
            layout="constrained",
            figsize=(3 * len(self.template_names), 3),
        )
        axes[0].set_ylim(ylim)

        if coefficients is None:
            self.logger.info("Using regression results to plot.")
            coefficients = self.coefficients

        centers = (self.edges[:, :-1] + self.edges[:, 1:]) / 2
        data_binned = jax.vmap(
            partial(
                _my_bincount,
                accumutalor=jnp.zeros(self.nbins + 2, dtype=float),
                weights=self.data * self.weight_model(coefficients),
            )
        )(self.data_digitized)
        chi2_arr = (
            self.normalization * data_binned / self.randoms_binned - 1
        ) ** 2 / self.error**2
        chi2_arr_noweights = (
            self.normalization * self.data_binned_noweights / self.randoms_binned - 1
        ) ** 2 / self.error**2

        for index, (coefficient_name, ax) in enumerate(
            zip(self.template_names, axes, strict=True)
        ):
            partial_chi2_noweights = (chi2_arr_noweights[index]).sum()
            ax.errorbar(
                centers[index],
                self.normalization
                * self.data_binned_noweights[index]
                / self.randoms_binned[index],
                self.error[index],
                fmt=".",
                label=f"χ² = {partial_chi2_noweights:.2f}/{self.nbins} = {partial_chi2_noweights / self.nbins:.2f}",
            )

            partial_chi2 = (chi2_arr[index]).sum()
            ax.errorbar(
                centers[index],
                self.normalization * data_binned[index] / self.randoms_binned[index],
                self.error[index],
                fmt=".",
                label=f"χ² = {partial_chi2:.2f}/{self.nbins} = {partial_chi2 / self.nbins:.2f}",
            )

            hist_data_syst, bins = jnp.histogram(
                self.template_values_data[index], bins=nbinsh
            )
            y = (
                hist_data_syst / hist_data_syst.max() * 0.3 * (ylim[1] - ylim[0])
                + ylim[0]
            )
            ax.stairs(y, bins)

            ax.set_xlabel(coefficient_name)
            ax.legend()
        fig.supylabel("Density fluctuations")
        fig.suptitle(title)

        return fig, axes


def _read_systematic_templates_stacked_alt(ra, dec, sys_tab, use_maps, nside, nest):
    pixnums = hp.ang2pix(nside, np.radians(-dec + 90), np.radians(ra), nest=nest)
    systematics_values = sys_tab[use_maps][pixnums]
    return np.vstack(
        [systematics_values[column_name] for column_name in systematics_values.colnames]
    )


def read_catalog(file_path: str, columns: list[str] | None = None) -> np.ndarray:
    """
    Read a catalog from a FITS file and optionally select specific columns.

    Parameters
    ----------
    file_path : str
        Path to the FITS file containing the catalog data.
    columns : list[str] | None, optional
        List of column names to select from the catalog. If None, all columns are returned.

    Returns
    -------
    np.ndarray
        Numpy structured array containing the catalog data. If `columns` is specified, only the selected columns are returned; otherwise, all columns are returned.

    Notes
    -----
    This function uses fitsio under the hood, but avoids asking it to load certain columns only due to performance issues.
    """
    whole_catalog = fitsio.read(file_path)  # load ALL columns
    if columns:
        return whole_catalog[columns]
    else:
        return whole_catalog


def columns_for_weight_scheme(
    weight_scheme: str | None,
    redshift_colname: str,
    tracer: str,
) -> tuple[set, set]:
    """
    Return the columns needed for the data and the randoms given the weight computation of ``lklmini.produce_imweights.produce_imweights``, as a tuple of sets.

    Parameters
    ----------
    weight_scheme : str or None
        Weighting scheme used in the linear regression, see documentation of ``lklmini.produce_imweights.produce_imweights``. If set to None, corresponds to the clustering catalog option.
    redshift_colname : str
        Name of the column containing redshift information ("Z" or "Z_not4clus").
    tracer : str
        Name of the tracer. Useful for the full catalog (columns needed in ``goodz_infull``).

    Returns
    -------
    tuple[set, set]
        Two sets listing the required columns for the data and randoms respectively.

    Raises
    ------
    KeyError
        In the event of an unrecognized ``weight_scheme`` or ``tracer``.
    """
    data_colnames = {
        "RA",
        "DEC",
        "PHOTSYS",
        redshift_colname,
    }
    randoms_colnames = {
        "RA",
        "DEC",
        "PHOTSYS",
    }

    match tracer:
        case "ELG":
            data_colnames_full = {"o2c"}
        case "LRG" | "BGS":
            data_colnames_full = {"ZWARN", "DELTACHI2"}
        case "QSO":
            data_colnames_full = set()
        case _:
            raise KeyError("Unrecognized tracer %s", tracer)

    match weight_scheme:
        case None:  # clustering catalog
            return (
                data_colnames | {"WEIGHT", "WEIGHT_FKP", "WEIGHT_SYS", "WEIGHT_ZFAIL"},
                randoms_colnames
                | {
                    "WEIGHT",
                    "WEIGHT_FKP",
                    "WEIGHT_SYS",
                    "WEIGHT_ZFAIL",
                    redshift_colname,
                },
            )
        case "fracz":  # 1/FRACZ_TILELOCID based completeness weights
            return (
                data_colnames
                | data_colnames_full
                | {"FRACZ_TILELOCID", "FRAC_TLOBS_TILES"},
                randoms_colnames,
            )
        case "wt":  # whole weight column
            return (
                data_colnames | data_colnames_full | {"WEIGHT"},
                randoms_colnames | {"WEIGHT"},
            )
        case "wtfkp":  # whole weight column plus FKP
            return (
                data_colnames | data_colnames_full | {"WEIGHT", "WEIGHT_FKP"},
                randoms_colnames | {"WEIGHT", "WEIGHT_FKP"},
            )
        case "wt_comp":  # WEIGHT_COMP
            return (
                data_colnames | data_colnames_full | {"WEIGHT_COMP"},
                randoms_colnames | {"WEIGHT_COMP"},
            )
        case _:
            raise KeyError("Weight scheme %s is not recognized.", weight_scheme)


def make_fit_maps_dictionary(
    default: list[str],
    except_when: list[tuple[str, tuple[float, float], list[str]]],
    regions: list[str],
    redshift_range: list[tuple[float, float]],
) -> dict[str, dict[str, list[str]]]:
    """
    Create a dictionary of systematics maps to feed to :py:func:`~alaeboss.produce_imweights.produce_imweights`.

    In the simplest case, all maps are the same for all photometric regions and redshift ranges. If you wish to differ from this (for example to follow the DESI DR1 LRG procedure) then define the exceptions here.
    Be wary that regions must match to the internal ones in :py:func:`~alaeboss.produce_imweights.produce_imweights`, ie ["N", "S"] for most tracers and ["N", "SnotDES", "DES"] for QSOs. This function itself does not make any checks.

    Parameters
    ----------
    default : list[str]
        List of maps to default to.
    except_when : list[tuple[str, tuple[float, float], list[str]]]
        List of conditions and corresponding alternative maps, under the form [(region, (z_min, z_max), [map1, map2, ...]), ...]. Set to ``None`` to keep the default everywhere.
    regions : list[str]
        List of all photometry regions.
    redshift_range : list[tuple[float, float]]
        List of all redshift ranges.

    Returns
    -------
    dict[str, dict[str, list[str]]]
        Dictionary of the form ``{region:{(z_min, z_max):fit_maps}}``.

    Examples
    --------
    In a regular usecase with QSOs, where no deviation from default is expected (dummy map names for brevity):

    >>> make_fit_maps_dictionary(default=["map1", "map2"], except_when=None, regions=["N", "SnotDES", "DES"], redshift_range=[(0.8, 1.3), (1.3, 2.1), (2.1, 3.5)])
    {'N': {(0.8, 1.3): ['map1', 'map2'], (1.3, 2.1): ['map1', 'map2'], (2.1, 3.5): ['map1', 'map2']}, 'SnotDES': {(0.8, 1.3): ['map1', 'map2'], (1.3, 2.1): ['map1', 'map2'], (2.1, 3.5): ['map1', 'map2']}, 'DES': {(0.8, 1.3): ['map1', 'map2'], (1.3, 2.1): ['map1', 'map2'], (2.1, 3.5): ['map1', 'map2']}}

    With LRGs, using a different set of maps in the North (dummy map names for brevity):

    >>> make_fit_maps_dictionary(default=["map1", "map2"], except_when=[("N", (0.4, 0.6), ["map1"]), ("N", (0.6, 0.8), ["map1"]), ("N", (0.8, 1.3), ["map1"])], regions=["N", "S"], redshift_range=[(0.4, 0.6), (0.6, 0.8), (0.8, 1.1)])
    {'N': {(0.4, 0.6): ['map1'], (0.6, 0.8): ['map1'], (0.8, 1.1): ['map1', 'map2'], (0.8, 1.3): ['map1']}, 'S': {(0.4, 0.6): ['map1', 'map2'], (0.6, 0.8): ['map1', 'map2'], (0.8, 1.1): ['map1', 'map2']}}

    With LRGs, using a different set of maps for each redshift bin in the North:

    >>> make_fit_maps_dictionary(default=["map1", "map2"], except_when=[("N", (0.4, 0.6), ["map1"]), ("N", (0.6, 0.8), ["map2"]), ("N", (0.8, 1.3), ["map3"])], regions=["N", "S"], redshift_range=[(0.4, 0.6), (0.6, 0.8), (0.8, 1.1)])
    {'N': {(0.4, 0.6): ['map1'], (0.6, 0.8): ['map2'], (0.8, 1.1): ['map1', 'map2'], (0.8, 1.3): ['map3']}, 'S': {(0.4, 0.6): ['map1', 'map2'], (0.6, 0.8): ['map1', 'map2'], (0.8, 1.1): ['map1', 'map2']}}
    """
    fit_maps_dictionary = {
        region: dict.fromkeys(redshift_range, default) for region in regions
    }
    if except_when is not None:
        for region, zrange, alternative_maps in except_when:
            fit_maps_dictionary[region][zrange] = alternative_maps

    return fit_maps_dictionary


def produce_imweights(
    # Input and output control
    data_catalogs: np.ndarray,
    randoms_catalogs: np.ndarray,
    is_clustering_catalog: bool,
    tracer_type: str,
    redshift_range: list[(float, float)],
    templates_maps_path_S: str,  # noqa: N803
    templates_maps_path_N: str,  # noqa: N803
    fit_maps: list[str] | dict[str, dict[str, list[str]]],
    output_directory: str | None,
    output_catalog_path: str | None,
    weight_scheme: str | None,
    output_column_name: str | None = "WEIGHT_IMLIN",
    save_summary_plots: bool = True,
    # Regression-specific arguments
    nbins: int = 10,
    tail: float = 0.5,
    # Miscellaneous
    logger: logging.Logger | None = None,
    loglevel: str = "INFO",
    templates_maps_nside: int = 256,
    templates_maps_nested: bool = True,
):
    """
    Perform linear regression to compute imaging systematics weights for a given tracer type, data catalog, random catalogs, set of maps.

    This function reads in a data catalog and associated random catalogs, applies selection criteria, loads imaging systematics templates, and performs regression to estimate and assign imaging weights to the data. The regression is performed separately for different photometric regions and redshift bins. Optionally, summary plots can be saved.

    Parameters
    ----------
    data_catalogs : np.ndarray
        Concatenated input data catalogs.
    randoms_catalogs : np.ndarray
        Concatenated input randoms catalogs.
    tracer_type : str
        Type of tracer (e.g., 'LRG', 'ELG_LOP', 'QSO').
    redshift_range : list of tuple of float
        List of (z_min, z_max) tuples defining redshift bins for regression.
    templates_maps_path_S : str
        Path to the South region imaging systematics templates file.
    templates_maps_path_N : str
        Path to the North region imaging systematics templates file.
    templates_maps_nside : int
        Nside for the template maps that are being loaded from this path. Default is 256.
    templates_maps_nested : bool
        Whether template maps are in the nested scheme. Default is True.
    fit_maps : list[str] | dict[str, dict[str, list[str]]]
        List of template map names to use in the regression. If list of strings, will use the same maps for all regions and redshift ranges. Otherwise, can provide a dictionary of the form ``{region:{(z_min, z_max):fit_maps}}``. Function :py:func:`~alaeboss.produce_imweights.make_fit_maps_dictionary` can help produce these.
    output_directory : str or None
        Directory where output plots and parameter files will be saved. If None, nothing is saved.
    output_catalog_path : str or None
        Path to the catalog where the output weights should be written. If None, no output is written.
    weight_scheme : str or None
        Which weights to apply on the data and randoms (typically to account for uncompleteness when regressing). The corresponding columns need to be available in the catalog.
            * ``fracz``: 1/(``FRACZ_TILELOCID`` * ``FRAC_TLOBS_TILES``) for the data, 1 for the randoms.
            * ``wt``: ``WEIGHT`` column from the catalog for the data and for the randoms
            * ``wtfkp``: FKP weights, *ie* ``WEIGHT`` * ``WEIGHT_FKP`` for both data and randoms.
            * ``wt_comp``: ``WEIGHT_COMP`` column for the data, 1 for the randoms.

        In all previous cases, if ``WEIGHT_ZFAIL`` is available, the weights will be multiplied by it. The standard for full catalogs is "fracz".
        If you are using **clustering** catalogs, *ie* if ``is_clustering_catalog`` is set to True, ``weight_scheme`` should be set to None as it will not be used.
    output_column_name : str or None, optional
        Name of the output column to store computed weights in the data catalog (default is "WEIGHT_IMLIN"). If set to None, the original catalog will be left untouched.
    save_summary_plots : bool, optional
        Whether to save summary plots of the regression results (default is True).
    nbins : int, optional
        Number of bins to use in regression preparation (default is 10).
    tail : float, optional
        Fraction of data to cut as outliers from each tail (default is 0.5).
    logger : logging.Logger, optional
        Logger object for logging progress and information (default is None). This will not log anything if set to None.
    loglevel : str, optional
        Logging level for the regressor's logger (default is "INFO"). Will not affect ``logger``'s level.

    Returns
    -------
    numpy.ndarray
        The function modifies the input data catalog in place by adding or updating the output_column_name
        with computed imaging weights, and writes regression parameters and plots to the output directory.
        The computed imaging weights are returned as an array (same shape as the input data, set to 1.0 where weights weren't computed.)

    Notes
    -----
    Loading some columns only from FITS file during NERSC jobs can be very long for mysterious reasons. If you are experiencing huge catalog readtimes, this might be why.

    """
    time_start = time()

    logger = logger or logging.getLogger("dummy")

    logger.info("Doing linear regression for imaging systematics")

    # Test that input parameters are compatible
    if is_clustering_catalog:
        logger.debug("The input catalogs are clustering catalogs...")
        if weight_scheme is not None:
            raise ValueError(
                "Cannot choose a weight scheme when using clustering catalogs ; `weight_scheme` should be set to `None`."
            )
        redshift_colname = "Z"
    else:
        assert weight_scheme is not None, (
            "Must define a weight scheme when using full catalogs."
        )
        redshift_colname = "Z_not4clus"

    # define photometric regions
    photometric_regions = ["S", "N"]
    if tracer_type == "QSO":
        photometric_regions = ["DES", "SnotDES", "N"]

    # Check if fit_maps is a list of strings or a dictionary
    if isinstance(fit_maps, dict):
        pass
    elif isinstance(fit_maps, Sequence) and isinstance(fit_maps[0], str):
        # do the conversion to a dictionary now
        fit_maps = make_fit_maps_dictionary(
            default=fit_maps,
            except_when=None,
            regions=photometric_regions,
            redshift_range=redshift_range,
        )
    else:
        raise TypeError(
            "Argument `fit_maps` is neither a dictionary nor a list of strings."
        )
    # Define a union of all fit maps to make sure to prepare any fit map needed
    all_fit_maps = list(
        {
            map_name
            for fitmap_r in fit_maps.values()
            for fitmap_rz in fitmap_r.values()
            for map_name in fitmap_rz
        }
    )
    all_fit_maps.sort()  # Keep everything in alphabetical order
    logger.debug("All fit maps: %s", all_fit_maps)

    # define a fit_maps dictionary in terms of subsets of all_fit_maps
    fit_maps_masks = {}  # For each region and redshift range, indicate which fit maps are used in terms of a mask of ``all_fit_maps``
    for region, fit_maps_r in fit_maps.items():
        fit_maps_masks[region] = {}
        for zrange, fit_maps_rz in fit_maps_r.items():
            fit_maps_masks[region][zrange] = np.array(
                [(mapname in fit_maps_rz) for mapname in all_fit_maps], dtype=bool
            )
    logger.debug("Masks map dictionary: %s", fit_maps_masks)

    jax.config.update("jax_enable_x64", True)
    logger.info("Enabled 64-bit mode for JAX")

    # define which columns need to be present for the data and the randoms
    data_colnames, random_colnames = columns_for_weight_scheme(
        weight_scheme=weight_scheme,
        redshift_colname=redshift_colname,
        tracer=tracer_type[:3],
    )
    data_colnames = list(data_colnames)
    random_colnames = list(random_colnames)
    logger.debug("Columns necessary for the data: %s", data_colnames)
    logger.debug("Columns necessary for the randoms: %s", random_colnames)

    debv = common.get_debv()  # for later
    sky_g, sky_r, sky_z = common.get_skyres()
    if output_directory is not None:
        output_directory = Path(output_directory)

    # copy data/randoms to avoid modifying in place
    # and check that necessary columns are present
    logger.info("There are %i rows of data", data_catalogs.size)
    all_data = Table(data_catalogs[data_colnames], copy=True)
    logger.info("There are %i rows of data", randoms_catalogs.size)
    rands = randoms_catalogs[random_colnames].copy()

    # select good data that has been observed
    if not is_clustering_catalog:
        logger.info("Selecting good and observed data")
        data_selection = common.goodz_infull(tracer_type[:3], all_data) & (
            all_data["ZWARN"] != 999999
        )
    else:
        data_selection = np.full_like(all_data, fill_value=True, dtype=bool)
    dat = all_data[data_selection]

    # prepare array to receive computed weights
    weights_imlin = np.ones(len(dat), dtype=float)

    for region in photometric_regions:
        # Are we north or south? (irrespective of DES)
        northsouth = region
        if region in ["DES", "SnotDES"]:
            northsouth = "S"

        # get healpix maps for the systematics
        logger.info("Loading healpix templates for region %s", northsouth)
        match northsouth:
            case "S":
                sys_tab = Table.read(templates_maps_path_S)
            case "N":
                sys_tab = Table.read(templates_maps_path_N)
            case _:
                raise KeyError(
                    "Value %s is not valid as a template region (North or South, ie 'S' or 'N')",
                    northsouth,
                )

        cols = list(sys_tab.dtype.names)  # names of templates

        for col in cols:
            if "DEPTH" in col:
                bnd = col.split("_")[-1]
                sys_tab[col] *= 10 ** (-0.4 * common.ext_coeff[bnd] * sys_tab["EBV"])
        for ec in ["GR", "RZ"]:
            sys_tab["EBV_DIFF_" + ec] = debv["EBV_DIFF_" + ec]
        if "EBV_DIFF_MPF" in all_fit_maps:
            sys_tab["EBV_DIFF_MPF"] = sys_tab["EBV"] - sys_tab["EBV_MPF_Mean_FW15"]
        if "SKY_RES_G" in all_fit_maps:
            sys_tab["SKY_RES_G"] = sky_g[northsouth]
        if "SKY_RES_R" in all_fit_maps:
            sys_tab["SKY_RES_R"] = sky_r[northsouth]
        if "SKY_RES_Z" in all_fit_maps:
            sys_tab["SKY_RES_Z"] = sky_z[northsouth]

        logger.info("All maps available for regression: %s", all_fit_maps)

        # select randoms now and retrieve systematics
        logger.info("Masking randoms according to the photometry region %s", region)
        match region:
            case "N" | "S":
                region_mask_randoms = rands["PHOTSYS"] == region
                region_mask_data = dat["PHOTSYS"] == region
            case "DES":
                region_mask_randoms = common.select_regressis_DES(rands)
                region_mask_data = common.select_regressis_DES(dat)
            case "SnotDES":
                region_mask_randoms = (rands["PHOTSYS"] == "S") & (
                    ~common.select_regressis_DES(rands)
                )
                region_mask_data = (dat["PHOTSYS"] == "S") & (
                    ~common.select_regressis_DES(dat)
                )
            case _:
                logger.info("other regions not currently supported")
                raise NotImplementedError("Exiting due to critical error with region")
        region_randoms = rands[region_mask_randoms]

        logger.info("Reading template values for the randoms")
        randoms_templates_values = _read_systematic_templates_stacked_alt(
            ra=region_randoms["RA"],
            dec=region_randoms["DEC"],
            sys_tab=sys_tab,
            use_maps=all_fit_maps,
            nside=templates_maps_nside,
            nest=templates_maps_nested,
        )
        logger.info(
            "Preparation for region %s is done. Starting regressions per redshift slice.",
            region,
        )

        for z_range in redshift_range:
            logger.info(
                "Getting weights for region %s and redshift bin %f < z < %f",
                region,
                z_range[0],
                z_range[1],
            )
            local_fit_maps = fit_maps[region][z_range]
            local_fit_maps.sort()  # force same order as all_fit_maps
            local_fit_maps_mask = fit_maps_masks[region][z_range]
            logger.info("Fit maps for this region are %s", local_fit_maps)
            t1 = time()
            # select data
            logger.info("Selecting data and loading template values")
            selection_data = (
                region_mask_data
                & (dat[redshift_colname] > z_range[0])
                & (dat[redshift_colname] < z_range[1])
            )
            selected_data = dat[selection_data]

            # if using clustering catalogs, select randoms further
            if is_clustering_catalog:
                logger.info("Clustering catalogs : selecting randoms")
                selection_randoms = (region_randoms[redshift_colname] > z_range[0]) & (
                    region_randoms[redshift_colname] < z_range[1]
                )
                selected_randoms = region_randoms[selection_randoms]
                selected_randoms_templates_values = randoms_templates_values[
                    local_fit_maps_mask, :
                ][:, selection_randoms]
            else:
                selected_randoms = region_randoms
                selected_randoms_templates_values = randoms_templates_values[
                    local_fit_maps_mask, :
                ]

            # get data imaging systematics
            data_templates_values = _read_systematic_templates_stacked_alt(
                ra=selected_data["RA"],
                dec=selected_data["DEC"],
                sys_tab=sys_tab,
                use_maps=local_fit_maps,
                nside=templates_maps_nside,
                nest=templates_maps_nested,
            )

            data_we = jax.numpy.ones_like(selected_data[redshift_colname], dtype=float)
            rand_we = jax.numpy.ones_like(selected_randoms, dtype=float)

            # add weights
            datacols = list(selected_data.dtype.names)
            logger.info("Found columns %s", datacols)

            match weight_scheme:
                case None:
                    assert is_clustering_catalog, (
                        "Cannot set weight_scheme to None if the catalogs are not clustering catalogs!"
                    )
                    logger.info(
                        "Clustering catalogs: using WEIGHT * WEIGHT_FKP / WEIGHT_SYS"
                    )
                    data_we *= (
                        selected_data["WEIGHT"]
                        * selected_data["WEIGHT_FKP"]
                        / selected_data["WEIGHT_SYS"]
                        / selected_data["WEIGHT_ZFAIL"]
                    )  # will be re-multiplied by WEIGHT_ZFAIL later
                    rand_we *= (
                        selected_randoms["WEIGHT"]
                        * selected_randoms["WEIGHT_FKP"]
                        / selected_randoms["WEIGHT_SYS"]
                        / selected_randoms["WEIGHT_ZFAIL"]
                    )
                case "fracz":
                    logger.info("Using 1/FRACZ_TILELOCID based completeness weights")
                    data_we /= selected_data["FRACZ_TILELOCID"]
                    if "FRAC_TLOBS_TILES" in datacols:
                        logger.info("Using FRAC_TLOBS_TILES")
                        data_we /= selected_data["FRAC_TLOBS_TILES"]
                    else:
                        logger.info("no FRAC_TLOBS_TILES")
                case "wt":
                    logger.info("Using the WEIGHT column directly")
                    data_we *= selected_data["WEIGHT"]
                    rand_we *= selected_randoms["WEIGHT"]
                case "wtfkp":
                    logger.info("Using FKP weights and WEIGHT")
                    data_we *= selected_data["WEIGHT"] * selected_data["WEIGHT_FKP"]
                    rand_we *= (
                        selected_randoms["WEIGHT"] * selected_randoms["WEIGHT_FKP"]
                    )
                case "wt_comp":
                    logger.info("Using WEIGHT_COMP column directly")
                    data_we *= selected_data["WEIGHT_COMP"]
                    rand_we *= selected_randoms["WEIGHT_COMP"]
                case _:
                    logger.warning("Weight scheme %s is not recognized.", weight_scheme)

            if "WEIGHT_ZFAIL" in datacols:
                logger.info("Adding redshift failure weights to data weights")
                data_we *= selected_data["WEIGHT_ZFAIL"]
            else:
                logger.info("No redshift failure weights")

            logger.info("Starting regression...")
            regressor = LinearRegressor.from_stacked_templates(
                data_weights=data_we,
                random_weights=rand_we,
                template_values_data=data_templates_values,
                template_values_randoms=selected_randoms_templates_values,
                template_names=local_fit_maps,
                loglevel=loglevel,
            )
            regressor.cut_outliers(tail=tail)
            regressor.prepare(nbins=nbins)
            optimized_parameters = regressor.regress_minuit()

            logger.info("Regression done!")
            logger.info("Optimized parameters are %s", optimized_parameters)

            if output_directory is not None:
                output_directory.mkdir(parents=True, exist_ok=True)
                output_loc = (
                    output_directory
                    / f"{tracer_type}_{region}_{z_range[0]:.1f}_{z_range[1]:.1f}_linfitparam_jax.txt"
                )
                logger.info("Writing to %s", output_loc)
                with open(output_loc, "w") as fo:
                    for par_name, par_value in optimized_parameters.items():
                        fo.write(f"{par_name} {par_value}\n")

                if save_summary_plots:
                    figname = (
                        output_directory
                        / f"{tracer_type}_{region}_{z_range[0]:.1f}_{z_range[1]:.1f}_linimsysfit_jax.png"
                    )
                    logger.info("Saving figure to %s", figname)
                    fig, _ = regressor.plot_overdensity(ylim=[0.7, 1.3])
                    fig.savefig(figname)
            else:
                logger.info(
                    "`output_directory` set to None; not writing plots or parameters."
                )

            weights_imlin[selection_data] = regressor.export_weights()
            t2 = time()
            logger.info(
                "Done with region %s and redshift bin %f < z < %s, took %f seconds",
                region,
                z_range[0],
                z_range[1],
                t2 - t1,
            )

    time_end = time()
    logger.info(
        "All linear regressions are done, took %i seconds.",
        int(time_end - time_start),
    )

    if (output_column_name is not None) and (output_catalog_path is not None):
        all_data[output_column_name] = 1.0
        all_data[output_column_name][data_selection] = weights_imlin
        logger.info("Writing to disk at %s", output_catalog_path)
        common.write_LSS_scratchcp(
            dat, str(output_catalog_path), logger=logger
        )  # LSS logging cannot handle Path objects
        return all_data[output_column_name]
    else:
        logger.info("Not writing results to disk.")
        all_data_weights = np.ones_like(all_data[redshift_colname], dtype=float)
        all_data_weights[data_selection] = weights_imlin
        return all_data_weights
