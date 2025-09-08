"""
Taken from https://github.com/domichbt/alaeboss/tree/42b44429696e07ddb227e3a89d9f9a17dcafa770/alaeboss by fusing `linear_regressor.py` and `produce_imweights.py`
"""

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

import logging
from functools import partial
from pathlib import Path
from time import time

import fitsio
import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from jax.scipy.optimize import minimize as jaxminimize
from jax.typing import ArrayLike

import LSS.common_tools as common
import healpy as hp


@jax.jit
def my_bincount(idx, accumutalor, weights):
    return accumutalor.at[idx].add(weights)[1:-1]


class LinearRegressor:
    constant = "constant"  # Name of the constant contribution in the weights
    bin_margin = (
        1e-7  # Extra space at the beginning and end of edges to avoid null-sized bins
    )

    logger = logging.getLogger(__qualname__)
    logger.setLevel(logging.DEBUG)
    stream_handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s [%(name)s] %(levelname)s | %(message)s", "%Y-%m-%d %H:%M:%S"
    )
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    def __init__(
        self,
        data_weights: ArrayLike,
        random_weights: ArrayLike,
        templates: dict[str : (jnp.array, jnp.array)],
        loglevel: str = "INFO",
    ):
        """
        Initialize a LinearRegressor instance for imaging systematics, by giving information on the data and randoms' properties. The regression can then be run by calling the `cut_outliers`, `prepare` and `regress` methods in that order. Once the regression has succeeded, weights can then be obtained by calling the `weight_model` method on the `coefficients` attribute.

        :param data_weights: External, non-imaging weights for each object of the data catalog (FKP, completeness...)
        :type data_weights: ArrayLike
        :param random_weights: External, non-imaging weights for each object of the randoms catalog (FKP, completeness...)
        :type random_weights: ArrayLike
        :param templates: Dictionnary containing the imaging template values for each objects of the data and randoms catalogs, under the form of {template_name: (data_values, random_values)}
        :type templates: dict[str: (jnp.array, jnp.array)]
        :param loglevel: Minimal logging level, set by default to INFO. Additional output can be obtained from DEBUG.
        :type loglevel: str
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
                "`jax_enable_x64` is set to False; this means that computations will be capped to float32. The regression will most likely fail and at least return unexpected results. Are you sure you meant to do this? Otherwise use `jax.config.update('jax_enable_x64', True)"
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
        self._idx_to_template_name = {
            idx: name for idx, name in enumerate(self.template_names)
        }

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
        Initialize a LinearRegressor instance for imaging systematics, by giving information on the data and randoms' properties. The regression can then be run by calling the `cut_outliers`, `prepare` and `regress` methods in that order. Once the regression has succeeded, weights can then be obtained by calling the `weight_model` method on the `coefficients` attribute.

        :param data_weights: External, non-imaging weights for each object of the data catalog (FKP, completeness...)
        :type data_weights: ArrayLike
        :param random_weights: External, non-imaging weights for each object of the randoms catalog (FKP, completeness...)
        :type random_weights: ArrayLike
        :param template_values_data: Array containing the template values for the data, for each template (shape = (# of templates, length of data))
        :type template_values_data: ArrayLike
        :param template_values_randoms: Array containing the template values for the randoms, for each template (shape = (# of templates, length of randoms))
        :type template_values_randoms: ArrayLike
        :param template_names: List of the names of the templates, in the same order as in `template_values_data` and  `template_values_randoms`
        :type template_names: list[str]
        :param loglevel: Minimal logging level, set by default to INFO. Additional output can be obtained from DEBUG
        :type loglevel: str
        :return: LinearRegressor instance, as would have been initialized with templates written in dictionaries (default `__init__`).
        :rtype: LinearRegressor
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
        Save the current instance data to disk. No regression or bin related information is saved. After loading, do not repeat cut_outliers as they are already removed.

        :param filepath: Path to the file where data should be saved.
        :type filepath: str
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

        :param cls: Description
        :param filepath: Description
        :type filepath: str
        :param loglevel: Description
        :type loglevel: str
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
        """
        For each template, remove the data and random that possess values in the extremal tails of the distribution. In total, `tail`% of the values are removed for each template.

        :param tail: Percentage of the extremal values to remove for each template.
        :type tail: float

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
                my_bincount,
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
                my_bincount,
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
            (
                (
                    self.normalization
                    * self.data_binned_noweights
                    / self.randoms_binned
                    - 1
                )
                ** 2
                / self.error**2
            )
        )
        self.logger.debug("Computed initial, correction-less chi2")

    def weight_model(self, coefficients: ArrayLike):
        """
        Return the weights computed from the model for the data. The returned array has the same shape as ``.

        :param coefficients: Parameters of the regression to weight each template, as well as the constant offset at index 0.
        :type coefficients: ArrayLike

        """
        return 1 / (
            1
            + coefficients[0]
            + jnp.dot(coefficients[1:], self.template_values_data_normalized)
        )

    def chi2(self, coefficients: ArrayLike):
        """
        Return the total cost function, defined for each template as the square of the ratio (weighted data histogram) to (weighted random histogram) minus one.
        Coefficients are used for the weight computation for the data.
        Use Poisson error of the ratio of weighted data (using coefficients) and weighted randoms.

        :param coefficients: Parameters of the regression to weight each template, as well as the constant offset at index 0.
        :type coefficients: ArrayLike

        """
        data_binned = jax.vmap(
            partial(
                my_bincount,
                accumutalor=jnp.zeros(self.nbins + 2, dtype=float),
                weights=self.data * self.weight_model(coefficients),
            )
        )(self.data_digitized)
        # error_everywhere = self.normalization * \
        #     jnp.sqrt(data_binned / self.randoms_binned**2 + data_binned**2 / self.randoms_binned**3) # model dependent error
        # Compute the chisquare over actual imaging templates (not the constant, which is mostly useful for the weight model)
        return jnp.sum(
            (
                (self.normalization * data_binned / self.randoms_binned - 1) ** 2
                / self.error**2
            )
        )

    def regress(self, guess: ArrayLike | None = None) -> dict[str:float] | None:
        """
        Find optimal coefficients by minimizing the cost function. Can provide an initial guess `guess`, otherwise all ones will be used.

        :param guess: Initial guess for the regression coefficients, of shape (1 + number of templates).
        :type guess: ArrayLike | None
        :return: Dictionnary of the fit coefficient for each template. The offset is referred to as 'constant'. If the minimization fails, None is returned.
        :rtype: dict | None
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
            return dict(zip([self.constant] + self.template_names, res.x))

    def regress_minuit(self, guess: ArrayLike | None = None) -> dict[str:float] | None:
        """
        Find optimal coefficients by minimizing the cost function using iminuit. Can provide an initial guess `guess`, otherwise all ones will be used.

        :param guess: Initial guess for the regression coefficients, of shape (1 + number of templates).
        :type guess: ArrayLike | None
        :return: Dictionnary of the fit coefficient for each template. The offset is referred to as 'constant'. If the minimization fails, None is returned.
        :rtype: dict | None
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
            return dict(zip([self.constant] + self.template_names, list(m.values)))

    def export_weights(self):
        """
        Return the imaging systematics weights using the parameters computed from the regression. The weights are also computed for data with extremal values of the systematics that were ignored in the regression.
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

    def mask(self):
        if hasattr(self, "extremal_data_mask"):
            return self.good_values_data.at[self.good_values_data].set(
                self.extremal_data_mask
            )
        else:
            return self.good_values_data

    # Additional methods to match original `Syst` class
    def get_subsample(self, subdata_mask: ArrayLike):
        """
        Docstring for get_subsample

        :param self: Instance to copy
        :param subdata_mask: Mask to apply on the data and templates. Note that these should have had bad values and outliers cut out, so the shape may be different from initial data fed to the initializer. Correspondance between original and current shape can be obtained from `LinearRegressor.mask()`.
        :type subdata_mask: ArrayLike
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
        ylim: list[float, float] = [0.75, 1.25],
        nbinsh: int = 50,
        title: str = None,
    ):
        """
        Create a subplot for each template, and plot the normalized histogram values before and after the regression. The initial distribution of the template on the data is overlaid on the bottom of the plot.

        :param coefficients: parameters for the weight model; if None, will default to the regression result.
        :type coefficients: ArrayLike | None
        :param ylim: Limits for the y-axis
        :type ylim: list[float, float]
        :param nbinsh: Number of bins for the overlaid histograms
        :type nbinsh: int
        :param title: Optional title
        :type title: str
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
                my_bincount,
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

        for index, (coefficient_name, ax) in enumerate(zip(self.template_names, axes)):
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


def read_systematic_templates_stacked_alt(ra, dec, sys_tab, use_maps, nside, nest):
    pixnums = hp.ang2pix(nside, np.radians(-dec + 90), np.radians(ra), nest=nest)
    systematics_values = sys_tab[use_maps][pixnums]
    return np.vstack(
        [systematics_values[column_name] for column_name in systematics_values.colnames]
    )


def produce_imweights(
    # Input and output control
    data_catalog_path: str,
    random_catalogs_paths: list[str],
    tracer_type: str,
    redshift_range: list[(float, float)],
    templates_maps_path_S: str,
    templates_maps_path_N: str,
    fit_maps: list[str],
    output_directory: str,
    output_column_name: str = "WEIGHT_IMLIN",
    save_summary_plots: bool = True,
    # Regression-specific arguments
    nbins: int = 10,
    tail: float = 0.5,
    # Miscellaneous
    logger: logging.Logger = None,
    loglevel: str = "INFO",
    templates_maps_nside: int = 256,
    templates_maps_nested: bool = True,
):
    """
    Wrapper function to perform linear regression to compute imaging systematics weights for a given tracer type, data catalog, random catalogs, set of maps.
    This function reads in a data catalog and associated random catalogs, applies selection criteria, loads imaging systematics templates, and performs regression to estimate and assign imaging weights to the data. The regression is performed separately for different photometric regions and redshift bins. Optionally, summary plots can be saved.

    Parameters
    ----------
    data_catalog_path : str
        Path to the input data catalog FITS file.
    random_catalogs_paths : list[str]
        List of paths to random catalogs FITS files.
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
    fit_maps : list[str]
        List of template map names to use in the regression.
    output_directory : str
        Directory where output plots and parameter files will be saved.
    output_column_name : str, optional
        Name of the output column to store computed weights in the data catalog (default is "WEIGHT_IMLIN").
    save_summary_plots : bool, optional
        Whether to save summary plots of the regression results (default is True).
    nbins : int, optional
        Number of bins to use in regression preparation (default is 10).
    tail : float, optional
        Fraction of data to cut as outliers from each tail (default is 0.5).
    logger : logging.Logger, optional
        Logger object for logging progress and information (default is None). This will not log anything if set to None.
    loglevel : str, optional
        Logging level for the regressor's logger (default is "INFO"). Will not affect `logger`'s level.
    Returns
    -------
    None
        The function modifies the input data catalog in place by adding or updating the output_column_name
        with computed imaging weights, and writes regression parameters and plots to the output directory.
    Notes
    -----
    Loading some columns only from FITS file during NERSC jobs can be very long for mysterious reasons. If you are experiencing huge catalog readtimes, this might be why.

    Correspondence between arguments and variables in `mkCat_main.py`:
    * `templates_maps_path_N = f"{lssmapdirout}{tpstr}_mapprops_healpix_nested_nside{nside}_N.fits"`
    * `templates_maps_path_S = f"{lssmapdirout}{tpstr}_mapprops_healpix_nested_nside{nside}_S.fits"`
    * `redshift_range = zrl`
    * `tracer_type = type`
    * `output_directory = dirout`
    * `fit_maps`: the maps you want to regress against. Might be usemaps if not None, otherwise `mainp.fit_maps_allebv` (LRGs) or `mainp.fit_maps`
    * `output_directory = dirout`

    You can setup a basic logger to pass to the function as
    ```
    logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        stream_handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "%(asctime)s [%(name)s] %(levelname)s | %(message)s", "%Y-%m-%d %H:%M:%S"
        )
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)
    ```

    """
    logger = logger or logging.getLogger("dummy")

    logger.info("Doing linear regression for imaging systematics")

    jax.config.update("jax_enable_x64", True)
    logger.info("Enabled 64-bit mode for JAX")

    time_start = time()

    debv = common.get_debv()  # for later
    sky_g, sky_r, sky_z = common.get_skyres()
    output_directory = Path(output_directory)

    # read data catalogs
    logger.info("Reading data catalogs")
    all_data = Table(fitsio.read(data_catalog_path))
    # read randoms catalogs (note that since we are reading a subset of columns, this can take a lot on time from a job, no idea why)
    logger.info("Reading %i randoms catalogs", len(random_catalogs_paths))
    rands = np.concatenate(
        [
            fitsio.read(
                random_catalog_path,
                columns=["RA", "DEC", "PHOTSYS"],
            )
            for random_catalog_path in random_catalogs_paths
        ]
    )

    # select good data that has been observed
    logger.info("Selecting good and observed data")
    data_selection = common.goodz_infull(tracer_type[:3], all_data) & (
        all_data["ZWARN"] != 999999
    )
    dat = all_data[data_selection]

    # prepare array to receive computed weights
    weights_imlin = np.ones(len(dat), dtype=float)

    # define photometric regions
    photometric_regions = ["S", "N"]
    if tracer_type == "QSO":
        photometric_regions = ["DES", "SnotDES", "N"]

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
        if "EBV_DIFF_MPF" in fit_maps:
            sys_tab["EBV_DIFF_MPF"] = sys_tab["EBV"] - sys_tab["EBV_MPF_Mean_FW15"]
        if "SKY_RES_G" in fit_maps:
            sys_tab["SKY_RES_G"] = sky_g[northsouth]
        if "SKY_RES_R" in fit_maps:
            sys_tab["SKY_RES_R"] = sky_r[northsouth]
        if "SKY_RES_Z" in fit_maps:
            sys_tab["SKY_RES_Z"] = sky_z[northsouth]

        logger.info(f"Maps for regression: {fit_maps}")

        # select randoms now and retrieve systematics
        logger.info("Masking randoms according to the region")
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
        # rand_syst = densvar.read_systematic_maps_alt(region_randoms['RA'], region_randoms['DEC'], sys_tab, use_maps)
        logger.info("Reading template values for the randoms")
        randoms_templates_values = read_systematic_templates_stacked_alt(
            ra=region_randoms["RA"],
            dec=region_randoms["DEC"],
            sys_tab=sys_tab,
            use_maps=fit_maps,
            nside=templates_maps_nside,
            nest=templates_maps_nested,
        )
        logger.info(
            f"Preparation for region {region} is done. Starting regressions per redshift slide."
        )

        for z_range in redshift_range:
            logger.info(
                f"Getting weights for region {region} and redshift bin {z_range[0]} < z < {z_range[1]}"
            )
            t1 = time()
            # select data
            logger.info("Selecting data and loading template values")
            selection_data = (
                region_mask_data
                & (dat["Z_not4clus"] > z_range[0])
                & (dat["Z_not4clus"] < z_range[1])
            )
            selected_data = dat[selection_data]

            # don't select randoms further because we're not using the clustering catalogs anyways

            # get data imaging systematics
            data_templates_values = read_systematic_templates_stacked_alt(
                ra=selected_data["RA"],
                dec=selected_data["DEC"],
                sys_tab=sys_tab,
                use_maps=fit_maps,
                nside=templates_maps_nside,
                nest=templates_maps_nested,
            )

            # add weights
            datacols = list(selected_data.dtype.names)
            logger.info(f"Found columns {cols}")
            logger.info("Using 1/FRACZ_TILELOCID based completeness weights")
            wts = 1 / selected_data["FRACZ_TILELOCID"]
            if "FRAC_TLOBS_TILES" in datacols:
                logger.info("Using FRAC_TLOBS_TILES")
                wts *= 1 / selected_data["FRAC_TLOBS_TILES"]
            else:
                logger.info("no FRAC_TLOBS_TILES")
            if "WEIGHT_ZFAIL" in datacols:
                logger.info("Using redshift failure weights")
                wts *= selected_data["WEIGHT_ZFAIL"]
            else:
                logger.info("no redshift failure weights")

            data_we = jax.numpy.array(wts)
            rand_we = jax.numpy.ones_like(region_randoms, dtype=float)

            logger.info("Starting regression...")
            regressor = LinearRegressor.from_stacked_templates(
                data_weights=data_we,
                random_weights=rand_we,
                template_values_data=data_templates_values,
                template_values_randoms=randoms_templates_values,
                template_names=fit_maps,
                loglevel=loglevel,
            )
            regressor.cut_outliers(tail=tail)
            regressor.prepare(nbins=nbins)
            optimized_parameters = regressor.regress_minuit()

            logger.info("Regression done!")
            logger.info(f"Optimized parameters are {optimized_parameters}")

            output_loc = (
                output_directory
                / f"{tracer_type}_{region}_{z_range[0]:.1f}_{z_range[1]:.1f}_linfitparam_jax.txt"
            )
            logger.info(f"Writing to {output_loc}")
            with open(output_loc, "w") as fo:
                for par_name, par_value in optimized_parameters.items():
                    fo.write(str(par_name) + " " + str(par_value) + "\n")

            if save_summary_plots:
                figname = (
                    output_directory
                    / f"{tracer_type}_{region}_{z_range[0]:.1f}_{z_range[1]:.1f}_linimsysfit_jax.png"
                )
                logger.info(f"Saving figure to {figname}")
                fig, axes = regressor.plot_overdensity(ylim=[0.7, 1.3])
                fig.savefig(figname)

            weights_imlin[selection_data] = regressor.export_weights()
            t2 = time()
            logger.info(
                f"Done with region {region} and redshift bin {z_range[0]} < z < {z_range[1]}, took {t2 - t1} seconds"
            )

    time_end = time()
    logger.info(
        "All linear regressions are done, took %i seconds. Now writing to disk",
        int(time_end - time_start),
    )

    all_data[output_column_name] = 1.0
    all_data[output_column_name][data_selection] = weights_imlin
    common.write_LSS_scratchcp(
        dat, str(data_catalog_path), logger=logger
    )  # LSS logging cannot handle Path objects
