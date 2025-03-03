"""This file defines two classes: MSM and MSMCollection. MSM contains functions for
building and analysing Markov State Models contained by MSMCollection in such a way
that they are comparable, e.g. using the same clusters and metastable states"""

from os import path as _path, remove as _remove
import pickle as _pickle
import pyemma.plots as _plots
from pyemma.coordinates import cluster_kmeans as _kmeans, assign_to_centers as _assign_to_centers
from pyemma.msm import timescales_msm as _timescales, bayesian_markov_model as _bayesian_msm
from pyemma import config as _pyemma_config
import numpy as _np
from matplotlib.cm import ScalarMappable as _scalarmappable
from matplotlib.colors import Normalize as _colornorm
from matplotlib.pyplot import subplots as _subplots
from seaborn import light_palette as _light_palette, color_palette as _color_palette
from subprocess import run as _subprocess
from scipy.optimize import curve_fit as _curve_fit
from ammo.utils._utils import __parse_time as _parse_time


_pyemma_config.show_progress_bars = False


class MSMCollection:
    """A collection of MSMs, intended for easier comparison.
    """

    def __init__(self, file=None):
        """
        Create or load an MSM collection.

        Parameters
        ----------
        file : str
            file containing a saved MSM collection. If None, an empty collection will be created
        """
        self._MSMs = {}
        self.clusters = None

        if file is not None:
            self.__load(file)

    def __getitem__(self, title):
        return self._MSMs[title]

    def __fix_titles(self, titles):
        """Fix MSM titles for going through MSMs

        Parameters
        ----------
        titles : None, str, [str]
            titles to parse

        Returns
        -------
        correct_titles : [str]
            titles in the correct format
        """
        if titles is None:
            correct_titles = [key for key in self._MSMs.keys()]
        elif type(titles) == str:
            correct_titles = [titles]
        elif type(titles) == list:
            correct_titles = titles
        else:
            raise ValueError('Please give MSM titles as str or list')

        return correct_titles

    def __load(self, file):
        with open(file, 'rb') as fl:
            unpickled = _pickle.load(fl)
        
        for msm in unpickled.get_all().values():
            self.add_msm(msm)
        self.clusters = unpickled.clusters
        
        return None

    def get_all(self):
        """
        Get all MSMs as a dictionary
        """
        return self._MSMs

    def save(self, file):
        """
        Save MSM collection to a file

        Parameters
        ----------
        file : str
            file to save to

        Returns
        -------
        None
        """
        with open(file, 'wb') as file:
            _pickle.dump(self, file)

    def add_msm(self, msm):
        """Add an MSM to the MSMs dictionary.
        
        If locations are provided for new MSM, data will also be loaded.
        
        Parameters
        ----------
        msm :  allostery.msm.MSM, str
            an MSM object or title for a new MSM
        
        Returns
        -------
        None
        """
        if isinstance(msm, MSM):
            self._MSMs[msm.title] = msm
        elif isinstance(msm, str):
            self._MSMs[msm] = MSM(msm)
        else:
            raise ValueError('Please provide an MSM object or a name for a new MSM')

        return None

    def remove_msm(self, title):
        """Remove an MSM from the MSMs dictionary
        
        Parameters
        ----------
        title : str
            title of the MSM to remove
        
        Returns
        -------
        None
        """
        self._MSMs.pop(title)

        return None

    def load_data(self, titles, locations, file_names, trajectories=_np.arange(1, 101), frames=5000, timestep='10 ps', features=None, missing='ignore', reload=False):
        """Read the featurised trajectory data into the MSMs
        
        Parameters
        ----------
        titles : [str]
            a list of titles of MSMs
        
        locations : [[str]]
            a list of locations for each MSM to be loaded
        
        reload : bool
            whether to reload data if it is already loaded
        
        trajectories : numpy.array, [int]
            a list of trajectory numbers in the locations
        
        frames : int
            number of frames to load
        
        file_names : [str]
            names of files containing featurised trajectory data

        reload : bool
            load data again if data already exists
        
        Returns
        -------
        None
        """
        print('Loading data...')
        # check titles and locations match length
        titles = self.__fix_titles(titles)
        if len(titles) != len(locations):
            raise IndexError('Title and location list lengths must match')
        for title, location in zip(titles, locations):
            print(f'{title} MSM')
            if title not in self._MSMs.keys():
                # create MSM
                msm = MSM(title)
                msm.load_data(location, file_names, trajectories, frames, timestep, features, missing)
                # put MSM into dict
                self.add_msm(msm)
            elif reload:
                # reload existing MSM data
                self._MSMs[title].load_data(location, trajectories, frames, file_names)
            print('-' * 20)
        print('...done.')

        return None

    def get_all_data(self, titles=None):
        """Return all the available loaded data
        
        Parameters
        ----------
        titles : str, [str]
            a list of titles of MSMs whose data to get
        
        Returns
        -------
        all_data : numpy.array, (total_n_frames, n_features)
            all_data points
        """
        titles = self.__fix_titles(titles)

        all_data = _np.array(_np.vstack([_np.vstack(self._MSMs[key].data) for key in titles]))
        return all_data

    def cluster(self, titles=None, n_clusters=100, max_iter=50, centers=None):
        """Kmeans clustering of all specified trajectory data. Creates self.clusters.
        
        Parameters
        ----------
        titles : str, [str]
            titles of MSMs whose data will be used for clustering
        
        n_clusters : int
            number of microstates
        
        max_iter : int
            maximum iterations for kmeans clustering
        
        centers : numpy.array(n_clusters, n_features)
            predefined centers to use
        
        Returns
        -------
        None
        """
        titles = self.__fix_titles(titles)
        cluster_data = self.get_all_data(titles)

        self.clusters = _kmeans(cluster_data, k=n_clusters, max_iter=max_iter, clustercenters=centers, keep_data=True)

        return None

    def assign_to_clusters(self, titles=None, clusters=None):
        """Assign all MSM trajectory data to clusters
        
        Parameters
        ----------
        titles : str, [str]
            titles of MSMs whose data will be assigned to clusters
        
        clusters : numpy.array(n_clusters, n_features)
            custom cluster centers to assign to. Otherwise will be assigned to clusters.clustercenters
        
        Returns
        -------
        None
        """
        if titles is None:
            titles = [key for key in self._MSMs.keys()]
        if clusters is None:
            clusters = self.clusters.clustercenters

        for key in titles:
            self._MSMs[key].assign_to_clusters(clusters)

        return None

    def compute_its(self, titles=None, plot=True,
                    lags_to_try=(1, 5, 10, 25, 50, 100, 250, 500, 750, 1000, 1500, 2000, 2500, 3000), nits=10,
                    errors='bayes', time_units=None):
        """Compute Implied Timescales for each MSM
        
        Also plot each with the same yaxis scale if required.
        
        Parameters
        ----------
        titles : str, [str]
            titles of MSMs to compute ITS for
        
        plot : bool
            whether to plot the ITS
        
        lags_to_try : [int]
            lag times in trajectory data step units to try for the ITS
        
        nits : int
            number of implied timescales to be computed
         
        errors : str
            whether to compute statistical uncertainties

        plot : bool
            plot the ITS

        time_units : str
            units for plotting, e.g. 'ns' or 'us'. If None, the default "steps" will be used

        
        Returns
        -------
        fig : matplotlib.figure.Figure
            figure containing ITS plot axes
        
        ax : np.array of AxesSubplot
            plots of each ITS
        """
        print('Computing ITS...')
        if titles is None:
            titles = [key for key in self._MSMs.keys()]

        for key in titles:
            print(key)
            self._MSMs[key].compute_its(lags_to_try, nits, errors, plot=False)

        print('...done.')

        if plot:
            n = len(titles)
            fig, ax = _subplots(n, figsize=(6, 4.5 * n))
            if n == 1:
                ax = [ax]
            limits = []
            for i, key in enumerate(titles):
                _plots.plot_implied_timescales(self._MSMs[key].its, ax[i])
                limits.append(ax[i].get_ylim())
                ax[i].set_title(key)

            ylim = (_np.array(limits)[:, 0].min(), _np.array(limits)[:, 1].max())
            for i in range(n):
                ax[i].set_ylim(ylim)

            # fix axes units if required
            if time_units is not None:
                for i in range(n):
                    factor = _parse_time(self._MSMs[titles[i]].timestep, time_units, output_type='number')
                    ax[i].set_xticklabels([int(val) for val in ax[i].get_xticks()*factor])
                    ax[i].set_xlabel(f'lag time / {time_units}')
                    ax[i].set_yticklabels([val for val in ax[i].get_yticks()*factor])
                    ax[i].set_ylabel(f'timescale / {time_units}')
            fig.tight_layout()
            return fig, ax

    def build_msms(self, lag_time, titles=None, verbose=False):
        """Build each MSM with a specified lag time.

        If remove_disconnected=True, the connected sets for each MSM will be compared
        and only the cluster centers that are sampled by all systems will be used. New
        MSMs will be built using only those clusters.

        Parameters
        ----------
        lag_time : int, str
            MSM lag time in trajectory steps (if int) or in format "value unit", e.g. "10 ps" (if str)

        titles : str, [str]
            systems to build MSMs for

        verbose : bool
            verbose output

        Returns
        -------
        None
        """
        titles = self.__fix_titles(titles)

        print('Building MSMs...')
        for key in titles:
            self._MSMs[key].build_msm(lag_time)
        print('...done.\n')

        return None

    def plot_data(self, shape, titles=None, x=0, y=1, features='infer', cmap=None):
        """Plot the data for each MSM as a density plot.

        All plots will use the same x and y axis scale, as well as the same
        density scale.

        Parameters
        ----------
        shape : tuple
            plot layout in the form (columns, rows)

        titles : str, [str]
            MSMs to plot the data of

        x : int
            x dimension

        y : int
            y dimension

        features : [str]
            a  list of x and y feature names for axis titles. "infer" will use "self.features" to find names. "None" will set no names

        cmap : matplotlib.colors.Colormap
            specify a colormap to use

        Returns
        -------
        fig : matplotlib.figure.Figure
            figure containing plot axes

        ax : np.array of AxesSubplot
            each MSM data as a density plot

        cbar : matplotlib.Colorbar
            colorbar
        """
        # fix inputs
        if cmap is None:
            cmap = _light_palette("seagreen", as_cmap=True)
        titles = self.__fix_titles(titles)

        fig, ax = _subplots(shape[0], shape[1], figsize=(5 * shape[1], 5 * shape[0]))
        if type(ax)==_np.ndarray:
            ax = ax.reshape(shape)
        else:
            ax = _np.array([ax]).reshape(shape)

        # find the same density limits
        clims = []
        for key in titles:
            X, Y, Z = _plots.plots2d.get_histogram(_np.vstack(self._MSMs[key].data)[:, x],
                                                         _np.vstack(self._MSMs[key].data)[:, y])
            density = _plots.plots2d._to_density(Z)
            clims.append((density.min(), density.max()))
        clims = (_np.array(clims)[:, 0].min(), _np.array(clims)[:, 1].max())

        # find the same axes limits
        limits = [_np.hstack([_np.vstack(self._MSMs[key].data)[:, x] for key in titles]).min(),
                  _np.hstack([_np.vstack(self._MSMs[key].data)[:, x] for key in titles]).max(),
                  _np.hstack([_np.vstack(self._MSMs[key].data)[:, y] for key in titles]).min(),
                  _np.hstack([_np.vstack(self._MSMs[key].data)[:, y] for key in titles]).max()]

        # reshape titles and plot
        titles = _np.array(titles).reshape(shape)
        for row in range(shape[0]):
            for col in range(shape[1]):
                fig_curr, ax_curr, misc_curr = _plots.plot_density(
                    _np.vstack(self._MSMs[titles[row, col]].data)[:, x],
                    _np.vstack(self._MSMs[titles[row, col]].data)[:, y],
                    ax=ax[row, col], cbar=False, vmax=clims[1], cmap=cmap)
                ax_curr.set_xlim((limits[0], limits[1]))
                ax_curr.set_ylim((limits[2], limits[3]))
                ax_curr.set_title(titles[row, col])
                if features == 'infer':
                    features = [self._MSMs[titles[row,col]].features[x], self._MSMs[titles[row,col]].features[y]]
                elif features == None:
                    features = [None, None]
                ax_curr.set_xlabel(features[0])
                ax_curr.set_ylabel(features[1])

        # plot colourbar legend
        cbar_ax = fig.add_axes([1.01, 0.18, 0.02, 0.7])
        cbar = fig.colorbar(
            _scalarmappable(norm=_colornorm(clims[0], clims[1]), cmap=cmap), cax=cbar_ax,
            label='density')
        fig.tight_layout()

        return fig, ax, cbar

    def pcca_assignments(self, n_states, msm=None, titles=None, disconnected=None, overwrite=False):
        """Compute metastable state probabilities from MSM stationary distributions

        Parameters
        ----------
        n_states : int
            number of metastable states

        msm : [str]
            MSMs whose pcca analysis to use. If None, all available will be used

        titles : str, [str]
            MSMs for which to compute the probabilities
            
        disconnected : int, [int]
            which metastable state to assign disconnected sets to. If list, length must match "msm" (i.e. different disconnected assignments for each pcca MSM). If None, will be assigned to the nearest metastable cluster
            
        overwrite : bool
            whether to overwrite existing probabilities

        Returns
        -------
        None
        """
        titles = self.__fix_titles(titles)
        if msm is None:
            msm = [key for key in self._MSMs.keys()]
        elif isinstance(msm, str):
            msm = [msm]
        if isinstance(disconnected, int) or disconnected is None:
            disconnected = [disconnected]*len(msm)
        if len(msm) != len(disconnected):
            raise ValueError(f'The length of "disconnected" ({disconnected}) must match "msm"')

        for key, state in zip(msm, disconnected):
            if n_states not in self._MSMs[key].pcca.keys() or overwrite:
                centers = self._MSMs[key].run_pcca(n_states, state)

            # compute assignments
            print(f'Metastable state assignments based on {key} MSM, {n_states} states:')
            for title in titles:
                print(title)
                self._MSMs[title].assign_to_metastable(n_states, state, self._MSMs[key], overwrite=overwrite)
                print('-' * 20)
            print()

        return None
    
    def add_state_labels(self, n_states, labels, titles=None):
        """Add meaningful labels for the PCCA assigned metastable states, e.g. "active" and "inactive"

        Parameters
        ----------
        n_states : int
            number of metastable states (for mapping to the corresponding PCCA assignments)

        labels : [str]
            list of state names. Length has to match the number of states

        titles: str, [str]
            MSMs to add the labels for. If None, will be added for all
        """
        titles = self.__fix_titles(titles)

        for title in titles:
            self._MSMs[title].add_state_labels(n_states, labels)

        return None

    def plot_stationary_distribution(self, shape, titles=None, x=0, y=1, features='infer', cmap=None, color='orange', plot_centres=True):
        """Plot stationary distributions of each microstate as a contour plot
        
        Parameters
        ----------
        shape : tuple
            plot layout in the form (columns, rows)

        titles : str, [str]
            MSMs to plot the data of

        x : int
            x dimension

        y : int
            y dimension

        features : [str]
            a  list of x and y feature names for axis titles. "infer" will use "self.features" to find names. "None" will set no names

        cmap : matplotlib.colors.Colormap
            specify a colormap to use

        color : str
            cluster color

 	plot_centres: bool
  	    decide whether to plot cluster centres
        """
        # fix inputs
        titles = self.__fix_titles(titles)

        fig, ax = _subplots(shape[0], shape[1], figsize=(5 * shape[1], 5 * shape[0]))
        if type(ax)==_np.ndarray:
            ax = ax.reshape(shape)
        else:
            ax = _np.array([ax]).reshape(shape)

        # find denisty limits
        probs = []
        for title in titles:
            probs += self._MSMs[title].stationary_distribution.tolist()
        clims = (min(probs), max(probs))

        # find the same axes limits
        limits = [_np.hstack([_np.vstack(self._MSMs[key].data)[:, x] for key in titles]).min(),
                  _np.hstack([_np.vstack(self._MSMs[key].data)[:, x] for key in titles]).max(),
                  _np.hstack([_np.vstack(self._MSMs[key].data)[:, y] for key in titles]).min(),
                  _np.hstack([_np.vstack(self._MSMs[key].data)[:, y] for key in titles]).max()]

        titles = _np.array(titles).reshape(shape)
        for row in range(shape[0]):
            for col in range(shape[1]):
                fig_curr, ax_curr, misc_curr = _plots.plot_contour(
                    _np.vstack(self._MSMs[titles[row, col]].data)[:, x],
                    _np.vstack(self._MSMs[titles[row, col]].data)[:, y],
                    _np.array(self._MSMs[titles[row, col]].stationary_distribution)[_np.hstack(self._MSMs[titles[row, col]].dtrajs)],
                    ax=ax[row, col], cbar=False, vmax=clims[1], cmap=cmap, method='nearest', mask=True)
                ax_curr.set_xlim((limits[0], limits[1]))
                ax_curr.set_ylim((limits[2], limits[3]))
                ax_curr.set_title(titles[row, col])
                if plot_centres:
                    ax_curr.scatter(self._MSMs[title].cluster_centers[:,x], self._MSMs[title].cluster_centers[:,y], s=8, c=color)
                if features == 'infer':
                    features = [self._MSMs[titles[row,col]].features[x], self._MSMs[titles[row,col]].features[y]]
                elif features == None:
                    features = [None, None]
                ax_curr.set_xlabel(features[0])
                ax_curr.set_ylabel(features[1])

        # plot colourbar legend
        cbar_ax = fig.add_axes([1.01, 0.18, 0.02, 0.7])
        cbar = fig.colorbar(
            _scalarmappable(norm=_colornorm(clims[0], clims[1]), cmap=cmap), cax=cbar_ax,
            label='stationary probability')
        fig.tight_layout()

        return fig, ax, cbar


    def plot_clusters(self, cluster_sets, shape, titles=None, x=0, y=1, features='infer', cmap=None, same_clusters=True, colors=['magenta', 'orange', 'teal', 'red']):
        """Plot the data for each MSM as a density plot, with cluster centers on top

        Parameters
        ----------
        cluster_sets : [np.array]
            a list of cluster sets to plot, e.g. 2 cluster sets of 50 clusters in 3 dimensions will have shape (2,50,3)
        
        shape : tuple
            plot layout in the form (columns, rows)

        titles : str, [str]
            MSMs to plot the data of

        x : int
            x dimension

        y : int
            y dimension

        features : [str]
            a  list of x and y feature names for axis titles. "infer" will use "self.features" to find names. "None" will set no names

        cmap : matplotlib.colors.Colormap
            specify a colormap to use

        same_clusters : bool
            whether to plot the same cluster sets on all data plots. If "False", it will be assumed that a correct number of individual cluster sets has been provided.

        colors : [str]
            cluster colors

        Returns
        -------
        fig : matplotlib.figure.Figure
            figure containing plot axes

        ax : np.array of AxesSubplot
            each MSM data as a density plot

        cbar : matplotlib.Colorbar
            colorbar
        """
        titles = self.__fix_titles(titles)
        fig, ax, cbar = self.plot_data(shape, titles, x, y, features, cmap)

        # sort the cluster sets
        if same_clusters: #if only one cluster set given
            cluster_sets = [cluster_sets for i in range(len(titles))]
        
        for row in range(shape[0]):
            for col in range(shape[1]):
                cluster = row*shape[1] + col
                for i, centers in enumerate(cluster_sets[cluster]):
                    ax[row, col].scatter(centers[:, x], centers[:, y], c=colors[i], s=10)

        return fig, ax, cbar

    def get_pcca_clusters(self, n_states, titles=None):
        """Get all PCCA cluster centers

        Parameters
        ----------
        n_states : int, [int]
            number of metastable states in the PCCA assignment. If a list of states, must
            match the number of titles given

        titles : str, [str]
            MSMs whose cluster enters to get

        Returns
        -------
        all_cluster_centers : numpy.array, shape (n_titles, 2)
            all cluster centers
        """
        titles = self.__fix_titles(titles)
        if type(n_states) == int:
            n_states = [n_states for i in range(len(titles))]
        elif len(n_states) != len(titles):
            raise IndexError('Lengths of titles and n_states must match')

        all_cluster_centers = []

        for key, states in zip(titles, n_states):
            msm = self._MSMs[key]
            centers = []
            for i in range(states):
                if len(msm.pcca[states][i]) > 0:
                    centers.append(msm.cluster_centers[msm.pcca[states][i]])
            all_cluster_centers.append(centers)

        return all_cluster_centers

    def mfpt(self, n_states, msm=None, timestep=None, titles=None, verbose=True, overwrite=False):
        """Compute mean first passage times for each MSM, based on specified pcca metastable state assignment

        Parameters
        ----------
        n_states : int
            number of metastable states

        msm : [str]
            MSMs to base the metastable assignment on. If None, all available will be used

        timestep : str
            trajectory timestep in format "value unit", e.g. "10 ps". If None, "self.timestep" will be used

        titles : str, [str]
            MSMs to compute MFPTs for

        verbose : bool
            verbose output

        overwrite : bool
            overwrite already computer MFPTs

        Returns
        -------
        None
        """
        titles = self.__fix_titles(titles)
        if msm is None:
            msm = [key for key in self._MSMs.keys()]

        for key in msm:
            # compute mfpts
            print(f'MFPTs based on {key} MSM, {n_states} states:')
            for title in titles:
                print(title)
                self._MSMs[title].compute_mfpt(n_states, timestep, msm=self._MSMs[key], overwrite=overwrite, verbose=verbose)
                print('-' * 20)
            print()

    def mfpr(self, n_states, msm=None, timestep=None, titles=None, verbose=True, overwrite=False):
        """Compute mean first passage rates for each MSM, based on specified pcca metastable state assignment

        Parameters
        ----------
        n_states : int
            number of metastable states

        msm : [str]
            MSMs to base the metastable assignment on. If None, all available will be used

        timestep : str
            trajectory timestep in format "value unit", e.g. "10 ps". If None, "self.timestep" will be used

        titles : str, [str]
            MSMs to compute MFPTs for

        verbose : bool
            verbose output

        overwrite : bool
            overwrite already computer MFPTs

        Returns
        -------
        None
        """
        titles = self.__fix_titles(titles)
        if msm is None:
            msm = [key for key in self._MSMs.keys()]

        for key in msm:
            # compute mfprs
            print(f'MFPRs based on {key} MSM, {n_states} states:')
            for title in titles:
                print(title)
                self._MSMs[title].compute_mfpr(n_states, timestep, msm=self._MSMs[key], overwrite=overwrite, verbose=verbose)
                print('-' * 20)
            print()

    def compare_states_and_timescales(self, n_states, msm=None, timestep=None, titles=None,
                                      overwrite=False):
        """Compare state probability ratios with MFPR ratios. These should be about the same if the system dynamics
        is well described by the metastable state assignment used.

        Parameters
        ----------
        n_states : int
            number of metastable states

        msm : [str]
            MSMs to base the metastable assignment on. If None, all available will be used

        timestep : str
            trajectory timestep in format "value unit", e.g. "10 ps". If None, "self.timestep" will be used

        titles : str, [str]
            MSMs to compute MFPTs for

        overwrite : bool
            overwrite already computer MFPTs

        Returns
        -------
        None
        """
        titles = self.__fix_titles(titles)
        if msm is None:
            msm = [key for key in self._MSMs.keys()]

        for key in msm:
            print(f'Timescales based on {key} MSM, {n_states} states:')
            for title in titles:
                print(title)
                self._MSMs[title].compare_states_and_timescales(n_states, timestep, msm=self._MSMs[key],
                                                               overwrite=overwrite)
                print('-' * 20)
            print()
        
        return None

    def bootstrapping(self, n_states, msm=None, titles=None, lag_time=None, cluster_centers=None, min_iter=100, max_iter=100, tol=1, last=10, verbose=False, overwrite=False):
        """
        Compute bootstrapped probabilities until they have converged to a Gaussian distribution or until maximum number
        of iterations have been reached.
        
        Parameters
        ----------
        n_states : int
            number of metastable states

        msm : allostery.msm.MSM
            MSMs whose metastable state assignment to use. If None, all will be used
        
        titles : str, [str]
            MSMs to compute probabilities for. If None, all with be used

        lag_time : int, str
            MSM lag time in trajectory steps (if int) or in format "value unit", e.g. "10 ps" (if str). If None, lag time of existing pyemma MSM will be used.
        
        cluster_centers : [float], numpy.array
            cluster centers to assign data to. If None, msm own cluster centers will be used
        
        min_iter : int
            minimum number of iterations
        
        max_iter : int
            maximum number of iterations
        
        tol : int, float
            tolerance for deviation from the mean
        
        last : int
            number of last resampling iterations that need to be within tolerance of the mean
        
        verbose : bool
            if verbose, progress will be reported

        overwrite : bool
            whether to overwrite existing probabilities
        
        Returns
        -------
        probabilities : {title: [float]}
            all of the computed active state probabilities as a dictionary
        """
        titles = self.__fix_titles(titles)
        if msm is None:
            msm = [key for key in self._MSMs.keys()]
        elif isinstance(msm, str):
            msm = [msm]
        if cluster_centers is None:
            cluster_centers = self.clusters.clustercenters
        
        probabilities = {}
        
        for pcca in msm:
            print(f'Bootstrapped probabilities, based on {pcca} MSM, {n_states} states:')
            for key in titles:
                print(key)
                probabilities[key] = self._MSMs[key].bootstrapping(n_states, self._MSMs[pcca], lag_time, cluster_centers, min_iter, max_iter, tol, last, verbose)
                print('-'*30)
        
        return probabilities
    
    def plot_bootstrapping_violin(self, assignment, xaxis='titles', titles=None, states=None, shape=None, colors=None):
        """Plot bootstrapped probabilities as a violin plot
        
        Parameters
        ----------
        assignment: str
            metastable state assingment to use, e.g. "reference, 2 states"

        xaxis : str
            what variable will be shown on the xaxis. Allowed values are: "titles" or "states". The other value will be kept constant and have to be provided as variables

        titles : str, [str]
            msm titles. If None, all will be used. If "xaxis=state", titles has to be a str of one MSM title

        states : str, int, [str], [int]
            metastable states to plot. If None, all will be plotted. Can be type str if state labels have been set, otherwise are integers starting from 1

        shape : tuple
            plot layout in the form (rows, columns)

        colors : [str]
            colors for each of the violin plots

        Returns
        -------
        fig : matplotlib.figure.Figure
            figure containing plot axes

        ax : np.array of AxesSubplot
            bootstrapped probabilities as a violin plot             
        """
        # metastable assignment breakdown
        assignment_title = assignment.split(',')[0]
        n_states = int(assignment.split(',')[1].split()[0])

        # clean titles and states input
        titles = self.__fix_titles(titles)
        if states is None:
            if n_states in self._MSMs[assignment_title].state_labels:
                states = self._MSMs[assignment_title].state_labels[n_states]
            else:
                states = [i+1 for i in range(n_states)]
        if not (isinstance(states, list) or isinstance(states, _np.ndarray)):
            states = [states]
        if isinstance(states[0], str):
            state_labels = states
            state_idxs = [_np.where(self._MSMs[assignment_title].state_labels[n_states]==state_name)[0][0] for state_name in states]
        elif isinstance(states[0], int):
            if n_states in self._MSMs[assignment_title].state_labels:
                state_labels = [self._MSMs[assignment_title].state_labels[3][i-1] for i in states]
            else:
                state_labels = [f'state {i}' for i in states]
            state_idxs = [i-1 for i in states]
        else:
            raise ValueError('"states" has to be int, str, or list of int or str')

        # get data and set axis labels depending on "xaxis"
        data = []
        if xaxis == 'titles':
            xticklabels = titles
            ylabels = [f'{state.capitalize()} state probability/%' for state in state_labels]
            for idx in state_idxs:
                data.append([self._MSMs[title].bootstrapping_data[assignment]['probabilities'][:,idx] for title in titles])
        elif xaxis == 'states':
            xticklabels = state_labels
            ylabels = [f'State probabilities of {title} MSM/%' for title in titles]
            for title in titles:
                data.append([self._MSMs[title].bootstrapping_data[assignment]['probabilities'][:,idx] for idx in state_idxs])

        # get figure parameters
        max_cols = 3
        if shape is None:
            shape = (len(data)//max_cols+int(len(data)%max_cols>0), max_cols)
        width = max(5, (len(data[0])*1.5))*shape[1]
        height = shape[0]*5

        # plot violin
        fig, ax = _subplots(shape[0], shape[1], figsize=(width, height))
        ax = _np.hstack(ax)
        for i in range(len(ax)):
            # if no more data, make invisible
            if i >= len(data):
                ax[i].set_visible(False)
            else:
                violins = ax[i].violinplot(data[i])
                # change color
                if colors is None:
                    colors = _color_palette(None, len(data[i]))
                for violin, color in zip(violins['bodies'], colors):
                    violin.set_color(color)
                violins['cmaxes'].set_color('black')
                violins['cbars'].set_color('black')
                violins['cmins'].set_color('black')
                # set ticks and labels
                ax[i].set_xticks(_np.arange(1, len(xticklabels)+1))
                ax[i].set_xticklabels(xticklabels, size=14)
                ax[i].set_ylabel(ylabels[i], size=14)

        fig.tight_layout()

        return fig, ax
    

class MSM:
    def __init__(self, title):
        self.stationary_distribution = None
        self.title = title
        self.locations = None
        self.timestep = None
        self.features = None
        self.__traj_locations = []
        self.data = None
        self.cluster_centers = None
        self.dtrajs = None
        self.its = None
        self.msm = None
        self.pcca = {}
        self.state_labels = {}
        self.metastable_assignments = {}
        self.metastable_assignments_bootstrapped = {}
        self.bootstrapping_data = {}
        self.mfpt = {}
        self.mfpr = {}

    def load_data(self, locations, file_names, trajectories=_np.arange(1, 101), frames=5000, timestep='10 ps', features=None, missing='ignore'):
        """Read the featurised trajectory data into the MSM
        
        Parameters
        ----------      
        locations : [[str]]
            a list of directories containing MSM data in folders called snapshot_[idx]

        file_names : [str]
            names of files containing featurised trajectory data
        
        trajectories : numpy.array, [int]
            a list of snapshot indices in 'locations'
        
        frames : int
            number of frames to load

        timestep : str
            trajectory timestep in format "value unit", e.g. "10 ps"

        features : [str]
            features that are being loaded

        missing : str
           how to handle missing data. Allowed values are "error", "warn", "ignore"
        
        Returns
        -------
        None
        """
        # set loactions
        self.locations = locations

        # set timestep
        self.timestep = timestep

        # set features
        if features is not None:
            if len(features) == len(file_names):
                self.features = features
            else:
                raise ValueError(f'"features" ({features}) should match "file_names" ({file_names})')
        else:
            self.features = [None for _ in range(len(file_names))]

        # load featurised data
        self.data = []
        missing_idx = []
        for directory in locations:
            print(directory)
            for idx in trajectories:
                print('trajectory %3i/%i' % (idx, trajectories[-1]), end='\r')
                # check if required files are there
                files = [f'{directory}/snapshot_{idx}/{name}' for name in file_names]
                if all([_path.exists(file) for file in files]):
                    #add to list
                    self.__traj_locations.append(f'{directory}/snapshot_{idx}')
                    # create data array
                    trajectory_data = _np.zeros((frames, 1))
                    for file in files:
                        trajectory_data = _np.append(trajectory_data, _np.loadtxt(file)[:frames].reshape(frames, 1),
                                                    axis=1)
                    # drop zeros and add to data
                    trajectory_data = _np.delete(trajectory_data, 0, 1)
                    self.data.append(trajectory_data)
                elif missing == 'error':
                    raise IOError(f'Data missing in {directory}/snapshot_{idx}')
                elif missing == 'warn':
                    missing_idx.append(f'{directory}/snapshot_{idx}')
        print()

        if len(missing_idx) > 0:
            print('Missing data:')
            print('\n'.join(missing_idx))

    def plot_data(self, x=0, y=1, features='infer', cmap=None):
        """Plot the data the MSM as a 2D density plot.

        Parameters
        ----------
        x : int
            x dimension

        y : int
            y dimension

        features : [str]
            a  list of x and y feature names for axis titles. "infer" will use "self.features" to find names. "None" will set no names

        cmap : matplotlib.colors.Colormap
            specify a colormap to use

        Returns
        -------
        fig : matplotlib.figure.Figure
            figure containing plot axes

        ax : np.array of AxesSubplot
            MSM data as a density plot

        misc : dict
            mappable and cbar
        """
        if cmap is None:
            cmap = _light_palette("seagreen", as_cmap=True)

        fig, ax, misc = _plots.plot_density(_np.vstack(self.data)[:, x], _np.vstack(self.data)[:, y], cmap=cmap)
        
        # get axis names
        if features == 'infer':
            features = [self.features[x], self.features[y]]
        elif features == None:
            features = [None, None]
        ax.set_xlabel(features[0])
        ax.set_ylabel(features[1])
        ax.set_title(self.title)

        # fix density bar to be more readable
        ticks = misc['cbar'].get_ticks()
        factor = int(f'{ticks[-1]:e}'.split('e')[-1])
        new_ticks = [round(tick*10**(-1*factor), 2) for tick in ticks]
        misc['cbar'].set_ticks(ticks)
        misc['cbar'].set_ticklabels(new_ticks)
        misc['cbar'].set_label(r'density $\times 10^%i$'%factor)

        return fig, ax, misc

    def cluster(self, n_clusters=100, max_iter=50, centers=None):
        """Kmeans clustering of all specified trajectory data. Creates self.clusters.
        
        Parameters
        ---------- 
        n_clusters : int
            number of microstates
        
        max_iter : int
            maximum iterations for kmeans clustering
        
        centers : numpy.array(n_clusters, n_features)
            predefined centers to use
        
        Returns
        -------
        None
        """
        self.cluster_centers = _kmeans(_np.vstack(self.data), k=n_clusters, max_iter=max_iter, clustercenters=centers, keep_data=True).clustercenters

        return None

    def assign_to_clusters(self, cluster_centers=None):
        """Assign the MSM trajectory data to clusters
        
        Parameters
        ----------
        cluster_centers : numpy.array
            cluster centers to assign to
            
        Returns
        ------
        None
        """
        if cluster_centers is not None:
            self.cluster_centers = cluster_centers
        self.dtrajs = _assign_to_centers(data=self.data, centers=self.cluster_centers)

    def compute_its(self, lags_to_try=(1, 5, 10, 25, 50, 100, 250, 500, 750, 1000, 1500, 2000, 2500, 3000), nits=10,
                    errors='bayes', plot=True, time_units=None):
        """Compute implied timescales

        Parameters
        ----------
        lags_to_try : [int]
            lags to try, in steps

        nits : int
            number of iterations

        errors : str
            error type

        plot : bool
            plot the ITS

        time_units : str
            units for plotting, e.g. 'ns' or 'us'. If None, the default "steps" will be used

        Return
        ------
        fig : matplotlib.figure.Figure
            figure containing plot axes

        ax : np.array of AxesSubplot
            ITS plot
        """
        self.its = _timescales(self.dtrajs, lags=lags_to_try, nits=nits, errors=errors)

        if plot:
            fig, ax = _subplots(1, figsize=(7,5))
            _plots.plot_implied_timescales(self.its, ax=ax)

            # fix axes units if required
            if time_units is not None:
                factor = _parse_time(self.timestep, time_units, output_type='number')
                ax.set_xticklabels([int(val) for val in ax.get_xticks()*factor])
                ax.set_xlabel(f'lag time / {time_units}')
                ax.set_yticklabels([val for val in ax.get_yticks()*factor])
                ax.set_ylabel(f'timescale / {time_units}')
            return fig, ax
        
        return None

    def build_msm(self, lag_time):
        """Build an MSM
        
        Parameters
        ----------
        lag_time : int, str
            MSM lag time in trajectory steps (if int) or in format "value unit", e.g. "10 ps" (if str)
            
        Returns
        -------
        None
        """
        if isinstance(lag_time, str):
            traj_step = _parse_time(self.timestep, 'ps', output_type='number')
            msm_step = _parse_time(lag_time, 'ps', output_type='number')
            lag_time = msm_step//traj_step
        self.msm = _bayesian_msm(self.dtrajs, lag_time)

        # add zero probability for disconnected sets
        self.stationary_distribution = self.msm.stationary_distribution
        disconnected_sets = self.msm.connected_sets[1:]
        if len(disconnected_sets)>0:
            disconnected_sets = _np.hstack(disconnected_sets)
            #sort disconnected sets so that the probabilities are inserted correctly
            disconnected_sets.sort()
            for center in disconnected_sets:
                self.stationary_distribution = _np.insert(self.stationary_distribution, center, 0.0)
                
        #if trajectory data was never assigned to a cluster centre at the end of the list, e.g. center 99, it does not appear in disconnected sets
        #therefore need to append 0.0 probability until there are the same number of stationary distributions as cluster centers
        while len(self.stationary_distribution)!=len(self.cluster_centers):
            self.stationary_distribution = _np.append(self.stationary_distribution, 0.0)

    def run_pcca(self, n_states, disconnected=None):
        """Run PCCA on the MSM
        
        Parameters
        ----------
        n_states : int
            number of states to assign
            
        disconnected : int
            which metastable state to assign disconnected sets to. If None, will be assigned to the nearest metastable cluster
            
        Returns
        -------
        pcca_centers
            cluster centers for each metastable state
        """
        # run pcca
        self.msm.pcca(n_states)

        # remap metastable sets to original clusters
        remapped = []
        for i in range(len(self.msm.metastable_sets)):
            metastable = [self.msm.connected_sets[0][center] for center in self.msm.metastable_sets[i]]
            remapped.append(_np.array(metastable))

        # ensure center indices are integers
        for state in range(n_states):
            remapped[state] = _np.array([int(center) for center in remapped[state]])

        # sort metastable sets based on average center distance from origin
        distances = _np.array([])
        for state in range(n_states):
            distances = _np.append(distances, _np.linalg.norm(self.cluster_centers[list(remapped[state])], axis=1).mean())
        # replace distances of 0 (i.e. empty sets) with highest value + offset
        # to put them at the back rather than the front
        distances[distances==0.0] = distances.max() + 10
        reordered = []
        for state in distances.argsort():
            reordered.append(remapped[state])

        # assign the disconnected sets to the specified
        # metastable state
        # if not specified, find the closest metastable cluster
        # to each disconnected cluster center
        if disconnected is None:
            cluster_av = _np.array([self.cluster_centers[list(reordered[i])].mean(axis=0) for i in range(n_states)])
            for i in range(1, len(self.msm.connected_sets)):
                for center_idx in self.msm.connected_sets[i]:
                    center_xyz = _np.array([self.cluster_centers[center_idx]])
                    to_add = _np.linalg.norm(cluster_av-center_xyz, axis=1).argmin()
                    modified_state = _np.append(reordered[to_add], center_idx)
                    modified_state.sort()
                    reordered[to_add] = modified_state.astype(int)
        else:
            state_to_add = reordered[disconnected]
            for i in range(1, len(self.msm.connected_sets)):
                for center_idx in self.msm.connected_sets[i]:
                    state_to_add = _np.append(state_to_add, center_idx)
            state_to_add.sort()
            reordered[disconnected] = state_to_add.astype(int)    

        self.pcca[n_states] = reordered

        # using list(centers) allows to give empty metastable sets
        pcca_centers = [self.cluster_centers[list(centers)] for centers in reordered]
        
        return pcca_centers
    
    def add_state_labels(self, n_states, labels):
        """Add meaningful labels for the PCCA assigned metastable states, e.g. "active" and "inactive"

        Parameters
        ----------
        n_states : int
            number of metastable states (for mapping to the corresponding PCCA assignments)

        labels : [str]
            list of state names. Length has to match the number of states
        """
        # check that number of states matches length
        if len(labels) != n_states:
            raise ValueError(f'Length of "labels" has to match "n_states"')
        
        self.state_labels[n_states] = _np.array(labels)

        return None

    def assign_to_metastable(self, n_states, disconnected=0, msm=None, verbose=True, overwrite=False, metastable_sets=None,
                             weights=None):
        """Compute metastable state probabilities

        Parameters
        ----------
        n_states : int
            number of metastable states

        disconnected : int
            which metastable state to assign disconnected sets to

        msm : allostery.msm.MSM
            MSMs whose pcca assignment to use. If None, own will be used

        verbose : bool
            whether to report results

        overwrite : bool
            whether to overwrite existing probabilities
            
        metastable_sets : numpy.array
            manual metastable sets
            
        weights : numpy.array
            manual weights for each state stationary probability. Must match metastable sets in shape

        Returns
        -------
        None
        """
        # use own pcca if none specified
        if msm is None:
            msm = self
        title = f'{msm.title}, {n_states} states'

        # use msm pcca metastable sets and weights
        # but can specify different sets and weights for modifying
        # which clusters should be included
        if metastable_sets is None:
            if n_states not in msm.pcca:
                msm.run_pcca(n_states, disconnected)
            metastable_sets = msm.pcca[n_states]
        if weights is None:
            weights = [_np.ones(len(sets)) for sets in metastable_sets]

        # if already assigned
        if title in self.metastable_assignments and not overwrite:
            if verbose:
                for key, info in self.metastable_assignments[title].items():
                    print(f'MS {key} has {info[2]} counts and {info[0]}% probability (± {info[1]}%)')
            return None
 
        # find the probabilities
        self.metastable_assignments[title] = {}
        for i in range(len(metastable_sets)):
            if len(metastable_sets[i]) > 0:
                probability = round((self.stationary_distribution[metastable_sets[i]] * weights[i]).sum() * 100, 2)
                error = round(_np.std(self.stationary_distribution[metastable_sets[i]] * weights[i]) * 100, 2)
                counts = len(metastable_sets[i])
            else:
                probability = 0
                error = 0
                counts = 0
            self.metastable_assignments[f'{msm.title}, {n_states} states'][i + 1] = [probability, error, counts]
            if verbose:
                print(f'MS {i + 1} has {counts} counts and {probability}% probability (± {error}%)')

    def plot_clusters(self, cluster_sets, x=0, y=1, features='infer', cmap=None, colors=['magenta', 'orange', 'teal', 'red']):
        """Plot cluster centers on top of data density plot
        
        Parameters
        ----------
        cluster_sets : [np.array]
            a list of cluster sets to plot, e.g. 2 cluster sets of 50 clusters in 3 dimensions will have shape (2,50,3)

        x : int
            x dimension

        y : int
            y dimension

        features : [str]
            a  list of x and y feature names for axis titles. "infer" will use "self.features" to find names. "None" will set no names

        cmap : matplotlib.colors.Colormap
            specify a colormap to use

        colors : [str]
            cluster colors

        Returns
        -------
        fig : matplotlib.figure.Figure
            figure containing plot axes

        ax : np.array of AxesSubplot
            MSM data as a density plot

        misc : dict
            mappable and cbar
        """
        fig, ax, misc = self.plot_data(x, y, features, cmap)
        for i, centers in enumerate(cluster_sets):
            ax.scatter(centers[:, x], centers[:, y], c=colors[i], s=8)

        return fig, ax, misc

    def plot_stationary_distribution(self, x=0, y=1, features='infer', cmap=None, color='orange'):
        """Plot stationary distributions of each microstate as a contour plot

        Parameters
        ----------
        x : int
            x dimension

        y : int
            y dimension

        features : [str]
            a  list of x and y feature names for axis titles. "infer" will use "self.features" to find names. "None" will set no names

        cmap : matplotlib.colors.Colormap
            specify a colormap to use

        color : str
            cluster color

        Returns
        -------
        fig : matplotlib.figure.Figure
            figure containing plot axes

        ax : np.array of AxesSubplot
            stationary distribution as a contour plot

        misc : dict
            mappable and cbar
        """
        z = _np.array(self.stationary_distribution)[_np.hstack(self.dtrajs)]
        fig, ax, misc = _plots.plot_contour(_np.vstack(self.data)[:,x], _np.vstack(self.data)[:,y], z, method='nearest', mask=True, cmap=cmap)
        ax.scatter(self.cluster_centers[:,x], self.cluster_centers[:,y], s=8, c=color)
        misc['cbar'].set_label('stationary probability')

        return fig, ax, misc

    def compute_mfpt(self, n_states, timestep=None, msm=None, verbose=True, overwrite=False):
        """Compute mean first passage times for the MSM, based on specified pcca metastable state assignment

        Parameters
        ----------
        n_states : int
            number of metastable states

        timestep : float
            trajectory timestep in format "value unit", e.g. "10 ps". If None, "self.timestep" will be used

        msm : allostery.msm.MSM
            MSM to base the metastable assignment on

        verbose : bool
            verbose output

        overwrite : bool
            overwrite already computer MFPTs

        Returns
        -------
        None
        """
        # get timestep
        if timestep is None:
            timestep = self.timestep
        step_time = _parse_time(timestep, 'us', 6, 'number')
        step_units = 'μs'

        # use own msm if none specified
        if msm is None:
            msm = self
        title = f'{msm.title}, {n_states} states'

        # if already computed
        if title in self.mfpt and not overwrite:
            if verbose:
                for key, info in self.mfpt[title].items():
                    print(
                        f'{key.split(",")[0]}->{key.split(",")[1]} transition: {round(info[0], 3)} {info[2]} (± {round(info[1], 3)} {info[2]})')
            return None

        # compute mfpt for the specified states
        self.mfpt[title] = {}
        # run pcca if missing
        if n_states not in msm.pcca:
            msm.run_pcca(n_states)
        # assigned sets
        metastable_sets = msm.pcca[n_states]

        # fix metastable sets to only include active sets
        for i, centers in enumerate(metastable_sets):
            new_centers = [_np.where(self.msm.active_set == center)[0][0] for center in centers if
                           center in self.msm.active_set]
            metastable_sets[i] = _np.array(new_centers)

        n_sets = len(metastable_sets)
        mfpt_sample = _np.zeros((n_sets, n_sets, self.msm.nsamples))
        # loop over sets
        for i in range(n_sets):
            for j in range(n_sets):
                if i != j:  # if a transition
                    mfpt_sample[i, j] = self.msm.sample_f('mfpt', metastable_sets[i], metastable_sets[j])
                    mfpt_sample[i, j] = mfpt_sample[i, j] * step_time  # get the right units
                    self.mfpt[title][f'{i + 1},{j + 1}'] = [mfpt_sample[i, j].mean(), _np.std(mfpt_sample[i, j]),
                                                            step_units]
                    if verbose:
                        print(
                            f'{i + 1}->{j + 1} transition: {round(mfpt_sample[i, j].mean(), 3)} {step_units} (± {round(_np.std(mfpt_sample[i, j]), 3)} {step_units})')

    def compute_mfpr(self, n_states, timestep=None, msm=None, verbose=True, overwrite=False):
        """Compute mean first passage rates for the MSM, based on specified pcca metastable state assignment

        Parameters
        ----------
        n_states : int
            number of metastable states

        timestep : float
            trajectory timestep in format "value unit", e.g. "10 ps". If None, "self.timestep" will be used

        msm : allostery.msm.MSM
            MSM to base the metastable assignment on

        verbose : bool
            verbose output

        overwrite : bool
            overwrite already computer MFPTs

        Returns
        -------
        None
        """
        # get timestep
        if timestep is None:
            timestep = self.timestep

        # use own msm if none specified
        if msm is None:
            msm = self
        title = f'{msm.title}, {n_states} states'

        # if already computed
        if title in self.mfpr and not overwrite:
            if verbose:
                for key, info in self.mfpr[title].items():
                    print(
                        f'{key.split(",")[0]}->{key.split(",")[1]} transition: {round(info[0], 3)}/{info[2]} (± {round(info[1], 3)}/{info[2]})')
            return None

        # compute mfpr from mfpt
        self.mfpr[title] = {}
        # compute mfpt with the same states and msm
        # if missing
        if title not in self.mfpt or overwrite:
            self.compute_mfpt(n_states, timestep, msm, False, overwrite)

        for key, info in self.mfpt[title].items():
            rate = 1 / info[0]
            error = rate * (info[1] / info[0])
            self.mfpr[title][key] = [rate, error, info[2]]
            if verbose:
                print(
                    f'{key.split(",")[0]}->{key.split(",")[1]} transition: {round(rate, 3)}/{info[2]} (± {round(error, 3)}/{info[2]})')

    def compare_states_and_timescales(self, n_states, timestep=None, msm=None, overwrite=False):
        """Compare state probability ratios with MFPR ratios. These should be about the same if the system dynamics
        is well described by the metastable state assignment used.

        Parameters
        ----------
        n_states : int
            number of metastable states

        timestep : float
            trajectory timestep in format "value unit", e.g. "10 ps". If None, "self.timestep" will be used

        msm : allostery.msm.MSM
            MSM to base the metastable assignment on

        overwrite : bool
            overwrite already computed  MFPTs

        Returns
        -------
        None
        """
        if msm is None:
            msm = self

        if f'{msm.title}, {n_states} states' not in self.metastable_assignments.keys():
            self.assign_to_metastable(n_states, msm, verbose=False, overwrite=overwrite)
        if f'{msm.title}, {n_states} states' not in self.mfpr.keys():
            self.compute_mfpr(n_states, timestep, msm, verbose=False, overwrite=overwrite)

        ratios = {}
        for i in range(1, n_states + 1):
            for j in range(i, n_states + 1):
                if i != j:
                    print(f'States {i} and {j}:')
                    states = [self.metastable_assignments[f'{msm.title}, {n_states} states'][i][0],
                              self.metastable_assignments[f'{msm.title}, {n_states} states'][j][0]]
                    transitions = [self.mfpr[f'{msm.title}, {n_states} states'][f'{i},{j}'][0],
                                   self.mfpr[f'{msm.title}, {n_states} states'][f'{j},{i}'][0]]
                    ratios[f'{i},{j}'] = (max(states) / min(states), max(transitions) / min(transitions))
                    print(f'State ratio: {round(max(states) / min(states), 2)}')
                    print(f'Transition timescale ratio: {round(max(transitions) / min(transitions), 2)}')
                    print()

        return ratios

    def sample_weighted_trajectories(self, n_frames, outputs=None, topology=None, distributions=None):
        """Sample trajectory data based on MSM distributions

        Parameters
        ----------
        n_frames : int
            number of frames to sample

        outputs : str, [str]
            output location. Number of outputs given must match number of distributions. If None, only samples will be returned

        topology : str
            topology file path. Can be None if "outputs" is None

        distributions : numpy.array
            distributions to base the sampling on. Must be same length as active set (i.e. removed disconnected states).
            If None, own stationary distribution will be used.

        Returns
        -------
        samples : numpy.array
            sampled trajectory and frame indices
        """
        if distributions is None:
            distributions = _np.array([self.msm.stationary_distribution])

        if type(outputs) == str:
            outputs = [outputs]

        if outputs is not None and len(outputs) != len(distributions):
            raise IndexError('Must provide the same number of outputs and distributions')

        samples = self.msm.sample_by_distributions(distributions, n_frames)

        # return here if no output
        if outputs is None:
            return samples
        
        # save if output provided
        for all_frames, output in zip(samples, outputs):
            file = '/tmp/cpptraj_sample_tmp.in'
            with open(file, 'w') as fl:
                fl.writelines(f'parm {_path.abspath(topology)}\n')
                for frame in all_frames:
                    fl.writelines(f'trajin {self.__traj_locations[frame[0]]}/production_dry.nc {frame[1]+1} {frame[1]+1}\n')
                fl.writelines(f'trajout {_path.abspath(output)}\n')
                fl.writelines('go\n')
                fl.writelines('quit\n')
            
            _subprocess(['cpptraj', '-i', file])
            _remove(file)

            # save which frames were sampled
            _np.savetxt(output.split('.')[0]+'.txt', all_frames, fmt='%i')

        return samples

    def __gauss(self, x, *p):
        A, mu, sigma = p
        return A*_np.exp(-(x-mu)**2/(2.*sigma**2))

    def __fit_gaus(self, data):
        """
        coeff : A, mu, sigma
        """
        data = _np.array(data)
        hist, bin_edges = _np.histogram(data)
        bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
        p0 = [10, data.mean(), 10]
        coeff, var_matrix = _curve_fit(self.__gauss, bin_centres, hist, p0=p0)

        return coeff

    def __build_bootstrapped_msm(self, lag, cluster_centers=None):
        """
        Build an msm with randomly resampled data and return state probability
        
        Parameters
        ----------
        lag : int
            msm lag time in steps        
        cluster_centers : [float], numpy.array
            cluster centers to assign data to. If None, msm own cluster centers will be used
        
        Returns
        -------
        stationary_distribution : [float]
            stationary distribution of bootstrapped msm

        traj_idxs : [int]
            indices of trajectories used for bootstrapped msm
        """
        #get resampled data
        traj_num = len(self.data)
        traj_idxs = _np.array([int(idx*traj_num) for idx in _np.random.rand(traj_num)])

	    # get new trajectories
	    # if different clusters provided, assign new dtrajs
        if cluster_centers is not None:
            new_data = [self.data[idx] for idx in traj_idxs]
            dtrajs = _assign_to_centers(new_data, cluster_centers)
	    # otherwise resample dtrajs directly
        else:
            dtrajs = [self.dtrajs[idx] for idx in traj_idxs]
            cluster_centers = self.cluster_centers
        #build msm
        bootstrap_msm = _bayesian_msm(dtrajs, lag)

        #get stationary distribution
        stationary_distribution = bootstrap_msm.stationary_distribution

        #add 0 probability for disconnected sets
        disconnected_sets = bootstrap_msm.connected_sets[1:]
        if len(disconnected_sets)>0:
            disconnected_sets = _np.hstack(disconnected_sets)
            disconnected_sets.sort()
            for center in disconnected_sets:
                stationary_distribution = _np.insert(stationary_distribution, center, 0.0)
        while len(stationary_distribution)!=len(cluster_centers):
            stationary_distribution = _np.append(stationary_distribution, 0.0)
            
        return stationary_distribution, traj_idxs

    def bootstrapping_convergence(self, state_probabilities, tol=1, last=10):
        """
        Check if bootstrapped probabilities have converged to a Gaussian distribution. This is done by
        checking whether the last n additions have changed the mean within some tolerance.
        
        Parameters
        ----------
        tol : int, float
            tolerance for devaition from the mean
        
        last : int
            number of last resampling iterations that need to be within tolerance of the mean
        
        Returns
        -------
        converged : bool
            whether bootstrapping has converged
        """
        # get the last mu values
        for i in range(state_probabilities.shape[1]):
            probabilities = state_probabilities[:,i]
            mu_values = []
            if len(probabilities)<last:#less than minimum needed to converge
                return False
            for i in range(last):
                last_idx = len(probabilities)-last+i+1
                try:
                    A, mu, sigma = self.__fit_gaus(probabilities[:last_idx])
                    mu_values.append(mu)
                except RuntimeError:#if cannot fit to a Gaussian curve, then not converged
                    return False
            mu_values = _np.array(mu_values)
            #check if converged
            av_value = mu_values.mean()
            diff = abs(mu_values-av_value)
            converged = all(diff<tol)
            if not converged:
                return False

        return converged

    def bootstrapping(self, n_states, msm=None, lag_time=None, cluster_centers=None, min_iter=10, max_iter=100, tol=1, last=10, verbose=False, overwrite=False):
        """
        Compute bootstrapped probabilities of a state until they have converged to a Gaussian distribution or until maximum number
        of iterations have been reached.
        
        Parameters
        ----------
        n_states : int
            number of metastable states

        msm : allostery.msm.MSM
            MSMs whose pcca assignment to use. If None, own will be used
        
        lag_time : int, str
            MSM lag time in trajectory steps (if int) or in format "value unit", e.g. "10 ps" (if str). If None, lag time of existing pyemma MSM will be used.
        
        cluster_centers : [float], numpy.array
            cluster centers to assign data to. If None, msm own cluster centers will be used
        
        min_iter : int
            minimum number of iterations
        
        max_iter : int
            maximum number of iterations
        
        tol : int, float
            tolerance for devaition from the mean
        
        last : int
            number of last resampling iterations that need to be within tolerance of the mean
        
        verbose : bool
            if verbose, progress will be reported

        overwrite : bool
            whether to overwrite existing probabilities
        
        Returns
        -------
        metastable_assignments_bootstraped: [float]
            mean computed probability and standard deviation
        """
        # prepare inputs
        # use own pcca if none specified
        if msm is None:
            msm = self

        # fix lag time
        if isinstance(lag_time, str):
            traj_step = _parse_time(self.timestep, 'ps', output_type='number')
            msm_step = _parse_time(lag_time, 'ps', output_type='number')
            lag_time = msm_step//traj_step
        elif lag_time is None:
            lag_time = self.msm.lagtime

        # check if bootstrapping is already done
        pcca = f'{msm.title}, {n_states} states'
        if pcca in self.bootstrapping_data and not overwrite:
            # set how many iterations already done
            # to add if more iterations are to be computed
            i = self.bootstrapping_data[pcca]['probabilities'].shape[0]+1
            converged = self.bootstrapping_convergence(self.bootstrapping_data[pcca]['probabilities'], tol, last)
        else:
            i = 1
            self.bootstrapping_data[pcca] = {'probabilities': _np.array([]), 'trajectories': []}
            converged = False
            
        # run bootstrapping
        while (not converged and i<=max_iter) or i<=min_iter:
            if verbose:
                print('%3i/%i'%(i,max_iter), end='\r')
            # build a bootstrapped msm
            try: 
                stationary_distribution, trajectories = self.__build_bootstrapped_msm(lag_time, cluster_centers)
            except Exception as e: # if msm stationary probabilities too low, an error is thrown - discard those
                continue
            # add results
            probability = _np.array([[round(stationary_distribution[list(state_clusters)].sum()*100, 2) for state_clusters in msm.pcca[n_states]]])
            if i == 1: # if first iteration
                self.bootstrapping_data[pcca]['probabilities'] = probability
            else:
                self.bootstrapping_data[pcca]['probabilities'] = _np.insert(self.bootstrapping_data[pcca]['probabilities'], len(self.bootstrapping_data[pcca]['probabilities']), probability, axis=0)
            self.bootstrapping_data[pcca]['trajectories'].append(trajectories)
            # check for convergence
            if i >= min_iter:
                converged = self.bootstrapping_convergence(self.bootstrapping_data[pcca]['probabilities'], tol, last)
            i+=1
            #except:
            #    None

        if verbose:
            print(' '*30, end='\r')
        if pcca not in self.metastable_assignments_bootstrapped or overwrite:
            self.metastable_assignments_bootstrapped[pcca] = {}
            for j in range(n_states):
                av_prob = round(_np.array(self.bootstrapping_data[pcca]['probabilities'][:,j]).mean(),2)
                sdev = round(_np.array(self.bootstrapping_data[pcca]['probabilities'][:,j]).std(), 2)
                self.metastable_assignments_bootstrapped[pcca][j+1] = [av_prob, sdev, i-1]
                print(f'State {j+1}: {av_prob}% ± {sdev}% ({i-1} iterations)')

        return self.metastable_assignments_bootstrapped[pcca]

    def plot_bootstrapping_violin(self, n_states, state_idx=1, msm=None, names=None, colors=None):
        """
        Plot bootstrapped probabilities as a violin plot

        Parameters
        ----------
        n_states : int
            number of metastable states in the metastable assignment

        state_idx : int, [int]
            index(es) in pcca assignment of state(s)

        msm : allostery.msm.MSM
            msm the metastable assignment was based on

        names : str, [str]
            state names, e.g. ["active", "inactive"]

        colors : [str]
            violin plot colors. Must be at least as long as "state_idx"

        Returns
        -------
        fig : matplotlib.figure.Figure
            figure containing plot axes

        ax : np.array of AxesSubplot
            bootstrapped probabilities as a violin plot
        """
        # prepare inputs
        # use own pcca if none specified
        if msm is None:
            msm = self
        # state indices as list
        if state_idx is None:
            state_idx = _np.arange(1, len(msm.pcca[n_states])+1)
        elif isinstance(state_idx, int):
            state_idx = [state_idx]
        # get state labels
        if names is None:
            names = [idx+1 for idx in state_idx]
        elif isinstance(names, str):
            names = [names]
        # get colors
        if colors is None:
            colors = _color_palette(None, n_states)

        # get data
        data = [self.bootstrapping_data[f'{msm.title}, {n_states} states']['probabilities'][:,idx-1] for idx in state_idx]

        # plot
        fig, ax = _subplots(1)
        violins = ax.violinplot(data)
        # change color
        for violin, color in zip(violins['bodies'], colors):
            violin.set_color(color)
            violin.set_alpha(0.6)
        violins['cmaxes'].set_color('black')
        violins['cbars'].set_color('black')
        violins['cmins'].set_color('black')
        # tidy ticks and labels
        ax.set_xticks(_np.arange(1, len(names)+1))
        ax.set_xticklabels(names)
        ax.set_xlabel('State')
        ax.set_ylabel('State probability/%')
        ax.set_ylim(0, None)

        return fig, ax
