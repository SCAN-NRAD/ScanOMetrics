"""
Module for implementation of developmental models
"""


import abc
import os
import numpy as np
from glob import glob1

from scanometrics import __version__
from scanometrics.utils import logging
from scipy.stats.mstats import mquantiles
from scipy.stats import norm
from scipy.stats import f as fstat
from scipy.stats import chi2
from numpy.polynomial import Polynomial as poly_model


class normative_model_template():

    def __init__(self, model_id, measured_metrics, metric_names, covariate_values, covariate_names,
                 dataset_id, proc_pipeline_id, cov2float, subject):
        """Constructor for the developmental model. Models are identified by unique combinations of model_name and
        dataset_id used to train it. The model folder is the folder to save/load fitted model parameters. Developmental
        models inherit ScanOMetrics_project attributes like measured_metrics, metric_names, covariates, subject dict,
        etc... that are used to fit the model and should be saved with it for reproducibility and plotting. The
        normative model has extra attributes like outliers, uncertainty, and stats with their respective methods to
        estimate them.

        :param dataset_id: unique ID of training dataset
        :type dataset_id: string
        """
        # Initialize attributes inherited from Scanometrics_project
        self.proc_pipeline_id = proc_pipeline_id  # ID of pipeline used for generate measured_metrics used to fit the model
        self.measured_metrics = measured_metrics.copy()
        self.metric_names = metric_names.copy()
        self.covariate_values = covariate_values.copy()
        self.covariate_names = covariate_names.copy()
        self.subject = subject.copy()
        self.cov2float = cov2float.copy()

        # Set own attributes
        self.model_id = model_id
        self.model_dataset_id = model_id + '_' + dataset_id   # Mixture of normative model (eg Polynomial) and training dataset ID (eg CHFIRST)

        # if subjects is None:  # Added at some point but probably not required: users should be able to load a trained model without downloading the training data.
        #     subjects = glob1(os.path.join(bids_folder, 'derivatives', 'scanometrics', 'sub-*'))
        # self.subjects = subjects
        self.parameter_fit = None  # estimated model parameters after fitting on self.measured_metrics (should be saved/loaded in/from model_folder/parameter_fit.npy)
        self.parameter_stat = None  # statistics of parameters trained on self.measured_metrics (should be saved/loaded in/from model_folder/parameter_stat.npy)
        self.outliers = {'orig': np.array([], dtype='bool'), 'norm': np.array([], dtype='bool')}
        self.uncertainty = {'orig': np.array([], dtype='float'), 'norm': np.array([], dtype='float')}
        self.stats = {}
        self.som_version = __version__

    def flag_outliers(self, k_IQR, variable='measured_metrics'):
        """
        Label subjects as outlier if morphological value :math:`x\\not\\in[q_{25}-k*IQR;q_{75}+k*IQR]`, where :math:`q_{25}` and :math:`q_{75}`
        are the 25th and 75th percentiles, IQR is the interquartile range, and k_IQR sets the threshold of how many IQRs are
        considered for labeling outliers. Subjects are compared to their age matching group [0.9*age,1.1*age]. The
        'metrics' matrix is a NxM matrix with N subjects and M metrics (eg self.measured_metrics or
        self.normativeModel.residuals). The quantiles are computed using numpy.quantile and the 'hazen' method, to
        obtain the same results as the Matlab 'quantile' function.

        :param k_IQR: factor of IQRs to be used as threshold for outlier detection.
        :type k_IQR: float > 0.
        """
        if k_IQR < 0:
            logging.ERROR('Parameter k_IQR in flag_outliers() is negative, should be >= 0.')
        if self.measured_metrics['orig'].shape == (0,):
            logging.ERROR('SOM project measured metrics is empty. Load metrics before running flag_outliers()')  # Adapt to residue test
        # Set residual or measured_metrics['orig'] as reference
        if variable == 'measured_metrics':
            metrics = self.measured_metrics  # Includes both orig and norm datasets
        elif variable == 'residuals':
            metrics = {}
            for k in self.fit_outputs.keys():
                metrics[k] = self.fit_outputs[k]['residuals']  # Should include both orig and norm datasets as keys
        else:
            logging.ERROR("flag_outliers(): %s variable not recognized (should be 'measured_metrics' or 'residuals')" %
                          variable)
        # Recover age vector
        age = self.covariate_values[:, self.covariate_names.index('age')]
        # Loop through orig and norm data
        self.stats = {}
        for metric_status in metrics.keys():
            # Initialize stat vectors
            self.stats[metric_status] = {'pout': np.zeros(metrics[metric_status].shape[1]),
                                         'fout': np.zeros(metrics[metric_status].shape[1]),
                                         'part': np.ones(metrics[metric_status].shape[1]),
                                         'odds': np.full(metrics[metric_status].shape[1], np.inf)}

            # If outliers has shape 0, meaning it's the first time flag_outliers is called, initialize array with measured_metrics.shape
            if self.outliers[metric_status].shape == (0,):
                self.outliers[metric_status] = np.zeros(metrics[metric_status].shape).astype('bool')
            # Copy current outliers : self.outliers is modified within the loop, which can lead to exclusion of subjects
            # depending on their order, and should be avoided
            outliers_preloop = self.outliers[metric_status].copy()
            # Loop through metrics m
            for m in range(metrics[metric_status].shape[1]):
                thr_lrg = np.full(metrics[metric_status].shape[0], np.nan)
                thr_sml = np.full(metrics[metric_status].shape[0], np.nan)
                # Loop through subjects s
                for s in range(metrics[metric_status].shape[0]):
                    if outliers_preloop[s, m]:  # skip current subject if already flagged as outlier
                        continue
                    # Find age matches
                    age_matches = ((age >= 0.9 * age[s]) & (age <= 1.1 * age[s]) & ~outliers_preloop[:, m])
                    age_matches[s] = False  # Remove current subject from matching array
                    if age_matches.sum() > 1:  # make sure there are at least two age matches for quantile computation
                        q_25, q_75 = np.quantile(metrics[metric_status][age_matches, m], [0.25, 0.75], method='hazen')
                        IQR = q_75 - q_25
                        thr_lrg[s] = q_75 + k_IQR * IQR
                        thr_sml[s] = q_25 - k_IQR * IQR
                        if metrics[metric_status][s, m] < thr_sml[s] or metrics[metric_status][s, m] > thr_lrg[s]:
                            self.outliers[metric_status][s, m] = True
                # Compute stats
                self.stats[metric_status]['fout'][m] = self.outliers[metric_status][:, m].sum() / metrics[metric_status].shape[0]  # compute outlier fraction
                if self.stats[metric_status]['fout'][m] > 0:
                    out_width = np.mean(np.abs(metrics[metric_status][self.outliers[metric_status][:, m], m]))
                    # We might run into issues with KIS variables according to Yujiang Wang (triggered complex valued array
                    # in octave implementation). For now this works:
                    self.stats[metric_status]['pout'][m] = norm.cdf(np.nanmean(thr_sml), loc=0, scale=out_width)\
                                   + 1 - norm.cdf(np.nanmean(thr_lrg), loc=0, scale=out_width)
                    if self.stats[metric_status]['pout'][m] > 0:
                        self.stats[metric_status]['part'][m] = self.stats[metric_status]['fout'][m] / self.stats[metric_status]['pout'][m]  # compute artifact probability
                        if self.stats[metric_status]['part'][m] > 0:
                            self.stats[metric_status]['odds'][m] = (1 - self.stats[metric_status]['part'][m]) / self.stats[metric_status]['part'][m]  # compute odds ratio
                    else:
                        self.stats[metric_status]['part'][m] = 1  # if self.pout is below or equal to zero,
                else:
                    self.stats[metric_status]['part'][m] = 0
                    self.stats[metric_status]['odds'][m] = np.inf

    def compute_uncertainty(self, min_repeat_subjs=2, frac=0.10):
        """
        Compute uncertainty over metrics in self.measured_metrics. Uncertainty is computed using repeated measures
        when available, otherwise approximates it with a given fraction of the mean across subjects (default frac=0.1).
        The uncertainty represents how much measured metrics are expected to deviate from their mean value, in the same
        units (eg mm for cortical thickness, mm^2 for surface area, etc...).
        Sets self.uncertainty to the computed value (vector with same number of values as columns in self.measured_metrics)
        Should be moved to normativeModel as it is currently copied from self. and not used outside context of normativeModel
        :param self: ScanOMmetrics project
        """
        # Check if there are repeats in subject list (within age*10%)
        rep_scan_indexes = []
        row_idx = 0
        # Get row indexes corresponding to subject
        logging.PRINT('Subjects with valid repeats:')
        for i_subj, ID in enumerate(self.subject.keys()):
            subj_rows = np.arange(row_idx, row_idx+len(self.subject[ID]))
            subj_ages = self.covariate_values[subj_rows, self.covariate_names.index('age')]
            sort_idx = np.argsort(subj_ages)
            i = 0  # start with youngest scan i
            while (i < len(sort_idx) - 1):  # keep i below n_scans -1, as checks for age @ i+1
                j = i + 1  # start comparison with scan j following the current reference i (j=i+1)
                valid_straw = [sort_idx[i]]  # Initialise straw with current scan (i)
                while (j<len(sort_idx) and subj_ages[sort_idx[j]]-subj_ages[sort_idx[i]]<0.1*subj_ages[sort_idx[i]]):  # if j was scanned before 0.1*age@j
                    valid_straw.append(sort_idx[j])
                    j += 1
                if len(valid_straw) > 1:
                    rep_scan_indexes.append(subj_rows[valid_straw])
                    logging.PRINT('\t-%s: %s' % (ID, [str(a) for a in self.covariate_values[subj_rows[valid_straw], self.covariate_names.index('age')]]))
                i += len(valid_straw)  # set reference scan to the one following the straw
            row_idx += len(self.subject[ID])
        if len(rep_scan_indexes) == 0:
            logging.PRINT('\t-[NONE]')
        # If repeats is None or empty, compute uncertainty as frac*mean(metric)
        self.uncertainty = {}
        for metric_status in self.measured_metrics.keys():
            if len(rep_scan_indexes) < min_repeat_subjs or len(self.outliers) == 0:
                self.uncertainty[metric_status] = frac*np.nanmean(self.measured_metrics[metric_status], 0)
            else:
                delta = np.zeros(self.measured_metrics[metric_status].shape[1])
                denom = np.zeros(self.measured_metrics[metric_status].shape[1])
                for idx in rep_scan_indexes:  # loop through rows corresponding to subjects with valid repeats
                    delta += (self.measured_metrics[metric_status].shape[0]-np.sum(self.outliers[metric_status][np.array(idx), :], 0)-1)*np.nanvar(self.measured_metrics[metric_status][np.array(idx), :], 0, ddof=0)  # repeats are usually few (2-3, rarely more than 10)
                    denom += self.measured_metrics[metric_status].shape[0]-np.sum(self.outliers[metric_status][np.array(idx), :], 0)
                self.uncertainty[metric_status] = np.sqrt(delta/(denom-1))


class Polynomial(normative_model_template):
    """
    Normative model based on polynomial fit. 'model_dataset_id' mixes normative
    model name with name of training dataset to keep track of combination used for training.
    """

    def __init__(self, measured_metrics, metric_names, covariate_values, covariate_names, dataset_id, proc_pipeline_id, cov2float, subject):
        """
        Class constructor. Initialize everything with empty arrays, matrices and lists.
        :param model_dataset_id: unique ID of dataset being processed.
        :type model_dataset_id: string
        """
        super().__init__('Polynomial', measured_metrics, metric_names, covariate_values, covariate_names,
                         dataset_id, proc_pipeline_id, cov2float, subject)
        self.age_vec          = None
        self.fit_outputs      = None

    def load_model_parameters(self):
        super().load_model_parameters()
        # Some model specific parameters

    def save_model_parameters(self):
        super().save_model_parameters()
        # Some model specific parameters

    def load_X(self):
        super().load_X()
        # Some model specific operation, although can probably be a general function without overriding...

    def save_X(self):
        super().save_X()
        # Some model specific operation, although can probably be a general function without overriding...

    def fit(self, flag_opt, global_deg_max=None, frac=20, alpha=0.01, N_cycl=10, n_uniform_bins=10):
        """Computes 'estimated_parameters' and 'residual' matrices for the model, according to 'measured_metrics' array.
        measured_metrics is the training dataset, given as a numpy array (a copy is created in self.measured_metrics to be saved/loaded). Outliers in
        the measured_metrics matrix still have values for computation of residuals.

        :param flag_opt: flag for optimizing maximum degree of polynomial. If set to False, the maximum degree is used.
                         If set to True, the SSE is computed for increasingly higher degrees, and stopped when there's
                         no significant improvement of the model.
        :type flag_opt: Bool
        """

        if 'age' not in self.covariate_names:
            logging.ERROR("Covariates in 'Polynomial' normative model requires at least an 'age' variable.")

        # Supposed to sample uniformly wrt age, but width of zero means all subjects are selected...
        # Setting uniform_selection to ones instead of calling uniform_subsample
        if N_cycl == 0:
            uniform_selection = np.ones((self.covariate_values.shape[0], 1), dtype='bool')
            N_cycl = 1
        else:
            uniform_selection = uniform_subsample_scheme(self.covariate_values[:, self.covariate_names.index('age')], N_cycl, n_uniform_bins)

        # Initialize empty arrays
        age      = self.covariate_values[:, self.covariate_names.index('age')].copy()
        self.age_vec = np.arange(np.floor(age.min()), np.ceil(age.max()) + 1)
        self.fit_outputs = {}
        for metric_status in self.measured_metrics.keys():
            if global_deg_max is None:
                deg_max = int(2 * np.floor(np.log(self.covariate_values.shape[0] / 10) + 1) - 1)  # according to https://www.frontiersin.org/articles/10.3389/fnins.2010.00058/full
            else:
                deg_max = global_deg_max
            self.fit_outputs[metric_status] = {'fit_deg': np.full((self.measured_metrics[metric_status].shape[1], N_cycl), np.nan),
                                               'fit_coeffs': np.zeros((deg_max + 1, self.measured_metrics[metric_status].shape[1], N_cycl)),
                                               'fit_good': np.zeros((self.measured_metrics[metric_status].shape[1], N_cycl)),
                                               'residuals': np.full(self.measured_metrics[metric_status].shape, np.nan),
                                               'fit_ave': np.full((len(self.age_vec), self.measured_metrics[metric_status].shape[1]), np.nan),
                                               'fit_dev': np.full((len(self.age_vec), self.measured_metrics[metric_status].shape[1]), np.nan),
                                               'est_dev': np.full((len(self.age_vec), self.measured_metrics[metric_status].shape[1]), np.nan),
                                               'fit_sml': np.full((len(self.age_vec), self.measured_metrics[metric_status].shape[1]), np.nan),
                                               'fit_lrg': np.full((len(self.age_vec), self.measured_metrics[metric_status].shape[1]), np.nan)
                                               }

            # Loop through the metrics
            for i in range(self.measured_metrics[metric_status].shape[1]):
                polyval_vec = np.zeros((len(self.age_vec), N_cycl))
                predict_vec = np.zeros((self.measured_metrics[metric_status].shape[0], N_cycl))
                for n in range(N_cycl):
                    idx_valid = np.argwhere(~self.outliers[metric_status][:, i] & ~np.isnan(self.measured_metrics[metric_status][:, i]) & uniform_selection[:, n])[:, 0]  # 0 indexing to get (N,) array TODO: consider changing to array of True/False
                    if len(idx_valid) > 2:
                        if global_deg_max is None:
                            deg_max = int(2 * np.floor(np.log(len(idx_valid)/10)+1)-1)  # according to https://www.frontiersin.org/articles/10.3389/fnins.2010.00058/full
                            if 'sym' in self.metric_names[i]:
                                deg_max = int(np.floor(deg_max))  # If assymetry index, take the floor integer
                            else:
                                deg_max = int(np.floor(deg_max) - (1 - np.floor(deg_max) % 2))  # If not assymetry, take lower odd number
                        else:
                            deg_max = global_deg_max
                        self.fit_outputs[metric_status]['fit_deg'][i, n] = deg_max
                        if flag_opt:
                            SSE = np.zeros(deg_max+1)  # to go from 1 to deg_max+1
                            try:
                                c = poly_model.fit(age[idx_valid], self.measured_metrics[metric_status][idx_valid, i], deg=0)
                            except np.linalg.LinAlgError:  # Tricky thing here is that it doesn't update sampling for previous metric_status if it converged
                                logging.WARNING('Convergence error, resampling uniform_selection, updating idx_valid and running fit again until it converges')
                                converged = False
                                while converged == False:
                                    logging.WARNING('Resampling...')
                                    uniform_selection[:, n] = uniform_subsample_scheme(self.covariate_values[:, self.covariate_names.index('age')], 1, n_uniform_bins)[:, 0]
                                    idx_valid = np.argwhere(~self.outliers[metric_status][:, i] & ~np.isnan(self.measured_metrics[metric_status][:, i]) & uniform_selection[:, n])[:, 0]
                                    while len(idx_valid) <= 2:
                                        idx_valid = np.argwhere(~self.outliers[metric_status][:, i] & ~np.isnan(self.measured_metrics[metric_status][:, i]) & uniform_selection[:, n])[:, 0]
                                    try:
                                        c = poly_model.fit(age[idx_valid], self.measured_metrics[metric_status][idx_valid, i], deg=0)
                                        converged = True
                                    except np.linalg.LinAlgError:
                                        continue
                            res = self.measured_metrics[metric_status][idx_valid, i] - c(age[idx_valid])
                            SSE[0] = res.dot(res.T)
                            for d in range(1, deg_max+1):  # goes from 1 to deg_max+1
                                try:
                                    c = poly_model.fit(age[idx_valid], self.measured_metrics[metric_status][idx_valid, i], deg=d)
                                except np.linalg.LinAlgError:  # Tricky thing here is that it doesn't update sampling for previous metric_status if it converged
                                    logging.WARNING('Convergence error, resampling uniform_selection, updating idx_valid and running fit again until it converges')
                                    converged = False
                                    while converged == False:
                                        logging.WARNING('Resampling...')
                                        uniform_selection[:, n] = uniform_subsample_scheme(self.covariate_values[:,self.covariate_names.index('age')], 1, n_uniform_bins)[:, 0]
                                        idx_valid = np.argwhere(~self.outliers[metric_status][:, i] & ~np.isnan(self.measured_metrics[metric_status][:,i]) & uniform_selection[:, n])[:, 0]
                                        while len(idx_valid) <= 2:
                                            idx_valid = np.argwhere(~self.outliers[metric_status][:, i] & ~np.isnan(self.measured_metrics[metric_status][:,i]) & uniform_selection[:, n])[:, 0]
                                        try:
                                            c = poly_model.fit(age[idx_valid],self.measured_metrics[metric_status][idx_valid, i], deg=d)
                                            converged = True
                                        except np.linalg.LinAlgError:
                                            continue
                                res = self.measured_metrics[metric_status][idx_valid, i] - c(age[idx_valid])
                                SSE[d] = res.dot(res.T)
                                # Stop if increase in d is not significant enough
                                df = len(idx_valid)-d
                                F = (SSE[d-1]-SSE[d]) / SSE[d] * df
                                pF = 1-fstat.cdf(F, 1, df)
                                if pF > alpha/d:
                                    self.fit_outputs[metric_status]['fit_deg'][i, n] = d-1
                                    break
                        if not np.isnan(self.fit_outputs[metric_status]['fit_deg'][i, n]):
                            d = int(self.fit_outputs[metric_status]['fit_deg'][i, n])
                            c = poly_model.fit(age[idx_valid], self.measured_metrics[metric_status][idx_valid, i], deg=d).convert()  # Gets coefficients for unscaled and unshifted polynomial basis, matching Matlab representation
                            self.fit_outputs[metric_status]['fit_coeffs'][:d+1, i, n] = c.coef.copy()
                            polyval_vec[:, n] = c(self.age_vec)
                            predict_vec[:, n] = c(age)
                        else:
                            logging.WARNING('Polynomial degree is NaN (using %s metric in %s measurements, i=%d' % (metric_status, self.metric_names[i], i))

                        # Chi-squared statistics to test whether the fit is good at all
                        # i.e. residual variance should not be significantly larger than variance of measurement uncertainty
                        T_fit = (len(idx_valid)-1) * np.var(self.measured_metrics[metric_status][idx_valid, i] - predict_vec[idx_valid, n], ddof=1) / self.uncertainty[metric_status][i]  # <- double check how uncertainty is computed
                        p_fit = 1-chi2.cdf(T_fit, len(idx_valid)-1)
                        if p_fit > alpha:
                            self.fit_outputs[metric_status]['fit_good'][i, n] = True
                idx_cycl = np.argwhere(self.fit_outputs[metric_status]['fit_good'][i, :])[:, 0]  # Indexes of n_cycl with good fit
                meas_is_good = ~np.isnan(self.measured_metrics[metric_status][:, i]) & ~self.outliers[metric_status][:, i]
                if len(idx_cycl) > 0:  # if there's at least one good fit, do the average over them
                    self.fit_outputs[metric_status]['fit_ave'][:, i] = polyval_vec[:, idx_cycl].mean(1)
                    if len(idx_cycl) > 1:
                        self.fit_outputs[metric_status]['fit_dev'][:, i] = polyval_vec[:, idx_cycl].std(1, ddof=1)
                    else:
                        self.fit_outputs[metric_status]['fit_dev'][:, i] = np.zeros(polyval_vec.shape[0])
                    for s in range(len(self.age_vec)):
                        idx_age = np.argwhere((age >= 0.9*self.age_vec[s]) & (age <= 1.1*self.age_vec[s]) & meas_is_good)[:, 0]
                        self.fit_outputs[metric_status]['est_dev'][s, i] = np.std(self.measured_metrics[metric_status][idx_age, i] - predict_vec[idx_age, :][:, idx_cycl].mean(1), ddof=1)
                else:
                    self.fit_outputs[metric_status]['fit_ave'][:, i] = np.nanmean(polyval_vec, 1)
                    if polyval_vec.shape[1] > 1:
                        self.fit_outputs[metric_status]['fit_dev'][:, i] = np.nanstd(polyval_vec, 1, ddof=1)
                    else:
                        self.fit_outputs[metric_status]['fit_dev'][:, i] = np.zeros(polyval_vec.shape[0])
                    for s in range(len(self.age_vec)):
                        idx_age = np.argwhere((age >= 0.9*self.age_vec[s]) & (age <= 1.1*self.age_vec[s]) & meas_is_good)[:, 0]
                        self.fit_outputs[metric_status]['est_dev'][s, i] = np.std(self.measured_metrics[metric_status][idx_age, i] - predict_vec[idx_age, :].mean(1), ddof=1)
                self.fit_outputs[metric_status]['fit_sml'][:, i] = self.fit_outputs[metric_status]['fit_ave'][:, i]-1.96*np.sqrt(self.fit_outputs[metric_status]['fit_dev'][:, i]**2+self.fit_outputs[metric_status]['est_dev'][:, i]**2)
                self.fit_outputs[metric_status]['fit_lrg'][:, i] = self.fit_outputs[metric_status]['fit_ave'][:, i]+1.96*np.sqrt(self.fit_outputs[metric_status]['fit_dev'][:, i]**2+self.fit_outputs[metric_status]['est_dev'][:, i]**2)
                if N_cycl > 0:  # original checks for size of predict_vec, should be equivalent
                    self.fit_outputs[metric_status]['residuals'][:, i] = self.measured_metrics[metric_status][:, i] - predict_vec.mean(1)

    def predict_values(self, new_covariate_values):
        # Loop through metrics to predict and compute residuals on
        predicted_values = {}
        for k in self.measured_metrics.keys():
            predicted_values[k] = np.full((new_covariate_values.shape[0], self.measured_metrics[k].shape[1]), np.nan)
            for i in range(self.measured_metrics[k].shape[1]):
                tmp = np.full((new_covariate_values.shape[0], self.fit_outputs[k]['fit_deg'].shape[1]), np.nan)
                for n in range(self.fit_outputs[k]['fit_deg'].shape[1]):
                    fit_deg = self.fit_outputs[k]['fit_deg'][i, n]
                    if ~np.isnan(fit_deg):
                        fit_deg = int(fit_deg)
                        c = poly_model.basis(fit_deg)
                        c.coef = self.fit_outputs[k]['fit_coeffs'][:fit_deg+1, i, n].copy()
                        for j in range(new_covariate_values.shape[0]):
                            tmp[j, n] = c(new_covariate_values[j, self.covariate_names.index('age')])
                if self.fit_outputs[k]['fit_good'][i, :].sum()>0:
                    predicted_values[k][:, i] = np.nanmean(tmp[:, self.fit_outputs[k]['fit_good'][i, :] > 0], axis=1)
                else:
                    predicted_values[k][:, i] = np.nanmean(tmp, axis=1)
        return predicted_values


def uniform_subsample_scheme_old(samples, n_cyl, width, seed=None):
    """Generates n_cyl subsampling schemes with approximate uniform distribution, by selecting subsamples with probability
    inversely proportional to the density.
    Original implementation from Octave's code. Possible shorter and pythonic implementation to be tested from here:
    https://stackoverflow.com/questions/66476638/downsampling-continuous-variable-to-uniform-distribution

    :param samples: sample from which subsamples should be taken.
    :param n_cyl: number of subsamples sets to generate
    :param width: width of broadening gaussian
    :param seed: seed for np.random.seed() for testing
    :return idx: np.array of size (len(sample),n_cyl), with 0 for excluded samples and 1 for included samples
    """
    if width > 0:
        density = np.zeros(len(samples))
        for s1 in range(len(samples)):
            for s2 in range(len(samples)):
                density[s1] += np.exp(-(samples[s2] - samples[s1]) ** 2 / (2 * width**2)) / (np.sqrt(2 * np.pi) * width)
        prob_sel = 1./density  # sum is approx equal to mean of sample ?
        # prob_sel /= prob_sel.sum()
        idx = (np.tile(prob_sel[:, None], (1, n_cyl)) > np.random.rand(len(samples), n_cyl)).astype('bool')
        # idx = np.random.choice(samples, len(samples), replace=True, p=prob_sel)
    else:
        idx = np.ones((len(samples), n_cyl)).astype('bool')
    return idx

def uniform_subsample_scheme(samples, n_cyl, n_bins=10):
    """Resample dataset to get an approximate uniform distribution (inspired from the following post on stackoverflow
    https://stackoverflow.com/questions/66476638/downsampling-continuous-variable-to-uniform-distribution). Taking the
    lowest bin count as sampling number can be too strict (i.e. it is likely that a bin has 1 or 2 subjects in it,
    which would lead to ~10 subjects selected for the analysis), but the method is quick and efficient.

    :param samples:
    :param n_cyl:
    :return:
    """
    counts, bins = np.histogram(samples, bins=n_bins)
    if np.amin(counts) == 0:
        logging.ERROR("Age distribution with %d bins leads to bin with zero subjects. The current uniform subsampling"
                      "scheme cannot handle such cases. Please reduce number of bins to get non-zero entries, or increase"
                      "your sample size." % n_bins)
    subsamples = np.full((len(samples), n_cyl), False)
    for cyl in range(n_cyl):
        for i in range(len(bins) - 1):
            if i == len(bins) - 2:
                # last bin is inclusive on both sides
                section = np.argwhere((samples >= bins[i]) & (samples <= bins[i + 1]))
            else:
                section = np.argwhere((samples >= bins[i]) & (samples < bins[i + 1]))
            sub_section = np.random.choice(section.reshape(-1), np.amin(counts), replace=False)
            subsamples[sub_section, cyl] = True
    return subsamples
