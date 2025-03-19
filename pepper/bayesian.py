import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.stats import lognorm
import emcee
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

class Bayesian:
    def __init__(self, y, comment_list = []):
        self.y = y
        self.set_comment_list(comment_list)
        # LOQ default settings
        self.LOQ_lower = -1   # (2.4 hours)
        self.LOQ_upper = 3 # 1000 days
        # prior default settings
        self.prior_mu_mean = 1
        self.prior_mu_std = 2
        self.prior_sigma_mean = 0.5
        self.prior_sigma_std = 0.5
        self.lower_limit_sigma = 0.05
        # EMCEE defaults
        self.nwalkers = 10
        self.iterations = 2000
        self.burn_in = 100
        self.save_backend = False
        # result variables
        self.sampler = None
        self.posterior_mu = None
        self.posterior_sigma = None
        self.posterior_mu_std= None
        # output options
        self.output_path = ''

    #### SET FUNCTIONS
    def set_comment_list(self, comment_list):
        if comment_list == [] or comment_list == None:
            self.comment_list = []
        else:
            assert type(comment_list) == list, "The comment list needs to be a list: {}".format(comment_list)
            assert len(comment_list) == len(self.y), "The size of the comment list provided does not match the size of y"
            self.comment_list = comment_list

    def set_prior_mu(self, mean, std):
        self.prior_mu_mean = mean
        self.prior_mu_std = std

    def set_prior_sigma(self, mean, std):
        self.prior_sigma_mean = mean
        self.prior_sigma_std = std

    def set_walkers(self, n_walkers):
        self.nwalkers = n_walkers

    def set_lower_limit_sigma(self, lower_limit_sigma):
        self.lower_limit_sigma = lower_limit_sigma

    def set_iterations(self, n_iterations):
        self.iterations = n_iterations

    def set_burn_in(self, discard_cutoff):
        self.burn_in = discard_cutoff

    def set_save_backend(self, save_backend):
        self.save_backend = save_backend

    def set_path_to_output_folder(self, output_path):
        self.output_path = output_path

    def set_LOQ_upper(self, LOQ):
        self.LOQ_upper = LOQ

    def set_LOQ_lower(self, LOQ):
        self.LOQ_lower = LOQ

    #### GET FUNCTIONS
    def get_walkers(self):
        return self.nwalkers

    def get_iterations(self):
        return self.iterations

    def get_burnin(self):
        return self.burn_in

    def get_LOQ_upper(self):
        return self.LOQ_upper

    def get_LOQ_lower(self):
        return self.LOQ_lower

    def get_prior_mu(self):
        return self.prior_mu_mean, self.prior_mu_std

    def get_prior_sigma(self):
        return self.prior_sigma_mean, self.prior_sigma_std

    def get_lower_limit_sigma(self):
        return self.lower_limit_sigma

    def get_sampler(self):
        if self.sampler is None:
            raise ValueError('Run sampling first before calling the sampler')
        return self.sampler

    def get_censored_values_only(self):
        censored_values = []
        for i, comment in enumerate(self.comment_list):
            if comment in ['<', '>']:
                censored_values.append(self.y[i])
            elif self.y[i] > self.LOQ_upper or self.y[i] < self.LOQ_lower:
                censored_values.append(self.y[i])
        return censored_values

    # Class functions
    def determine_LOQ(self):
        """
        Determines if the LOQ is upper or lower, and the value (if not default)
        :return: upper_LOQ , lower_LOQ
        """

        censored_values = self.get_censored_values_only()

        # Find upper LOQ
        upper_LOQ = np.nan
        # bigger than global LOQ
        if max(self.y) >= self.LOQ_upper:
            upper_LOQ = self.LOQ_upper
        # case if exactly 365 days
        elif max(self.y) == 2.562: # 365 days
            upper_LOQ = 2.562
            self.LOQ_upper = upper_LOQ
        # case if "bigger than" indication in comments
        elif '>' in self.comment_list:
            i = 0
            while i < len(self.y):
                if self.y[i] == min(censored_values) and self.comment_list[i] == '>':
                    self.LOQ_upper = self.y[i]
                    break
                i+=1

        # Find lower LOQ
        lower_LOQ = np.nan
        # smaller than global LOQ
        if min(self.y) <= self.LOQ_lower:
            lower_LOQ = self.LOQ_lower
        # case if exactly 1 day
        elif min(self.y) == 0: # 1 day
            lower_LOQ = 0
            self.LOQ_lower = 0
        # case if "smaller than" indication in comments
        elif '<' in self.comment_list:
            i = 0
            while i < len(self.y):
                if self.y[i] == max(censored_values) and self.comment_list[i] == '<':
                    self.LOQ_lower = self.y[i]
                    break
                i+=1
        return upper_LOQ, lower_LOQ

    def logLikelihood(self, theta, sigma):
        """
        Likelihood function (the probability of a dataset (mean, std) given the model parameters)
        Convert not censored observations into type ’numeric’
        :param theta: mean half-life value to be evaluated
        :param sigma: std half-life value to be evaluated
        :return: log_likelihood
        """
        upper_LOQ, lower_LOQ = self.determine_LOQ()

        n_censored_upper = 0
        n_censored_lower = 0
        y_not_cen = []

        if np.isnan(upper_LOQ) and np.isnan(lower_LOQ):
            y_not_cen = self.y
        else:
            for i in self.y:
                if np.isnan(upper_LOQ) and i >= upper_LOQ: # censor above threshold
                    n_censored_upper +=1
                if np.isnan(lower_LOQ) and i <= lower_LOQ: # censor below threshold
                    n_censored_lower += 1
                else: # do not censor
                    y_not_cen.append(i)

        LL_left_cen = 0
        LL_right_cen = 0
        LL_not_cen = 0

        ## likelihood for not censored observations
        if n_censored_lower > 0: # loglikelihood for left censored observations
            LL_left_cen = n_censored_lower * norm.logcdf(lower_LOQ, loc=theta, scale=sigma)  # cumulative distribution function CDF

        if n_censored_upper > 0: # loglikelihood for right censored observations
            LL_right_cen = n_censored_upper * norm.logsf(upper_LOQ, loc=theta, scale=sigma) # survival function (1-CDF)

        if len(y_not_cen) > 0: # loglikelihood for uncensored values
            LL_not_cen = sum(norm.logpdf(y_not_cen, loc=theta, scale=sigma)) # probability density function PDF

        return (LL_left_cen + LL_not_cen + LL_right_cen)

    def get_prior_probability_sigma(self, sigma):
        # convert mean and sd to logspace parameters, to see this formula check
        # https://en.wikipedia.org/wiki/Log-normal_distribution under Method of moments section
        temp = 1 + (self.prior_sigma_std / self.prior_sigma_mean) ** 2
        meanlog = self.prior_sigma_mean / np.sqrt(temp)
        sdlog = np.sqrt(np.log(temp))
        # calculate of logpdf of sigma
        norm_pdf_sigma = lognorm.logpdf(sigma, s=sdlog, loc=self.lower_limit_sigma,
                                        scale=meanlog)
        return norm_pdf_sigma

    def get_prior_probability_theta(self, theta):
        norm_pdf_theta = norm.logpdf(theta, loc=self.prior_mu_mean, scale=self.prior_mu_std)
        return norm_pdf_theta

    def logPrior(self, par):
         """
         Obtain prior loglikelihood of [theta, sigma]
         :param par: par = [theta,sigma]
         :return: loglikelihood
         """
         # calculate the mean and standard deviation in the log-space
         norm_pdf_mean = self.get_prior_probability_theta(par[0])
         norm_pdf_std = self.get_prior_probability_sigma(par[1])
         log_norm_pdf = [norm_pdf_mean, norm_pdf_std]
         return sum(log_norm_pdf)

    def logPosterior(self, par):
        """
        Obtain posterior loglikelihood
        :param par: [theta, sigma]
        :return: posterior loglikelihood
        """
        logpri = self.logPrior(par)
        if not np.isfinite(logpri):
            return -np.inf
        loglikelihood = self.logLikelihood(par[0], par[1])
        return logpri + loglikelihood

    def get_posterior_distribution(self):

        """
        Sample posterior distribution and get median of mean and std samples
        :return: posterior half-life mean and std
        """
        if self.posterior_mu:
            return self.posterior_mu, self.posterior_sigma

        # Sampler parameters
        ndim = 2 # number of dimensions (mean,std)
        p0 = abs(np.random.randn(self.nwalkers, ndim))  # only positive starting numbers (for std)
        # if log-prob data will be re-used later, it needs to be saved in backend file here
        backend = None
        if self.save_backend:
            backend = emcee.backends.HDFBackend(self.output_path + "backend.h5")
            backend.reset(self.get_walkers(), ndim)

        # Sample distribution
        print('Sampling in progress...')
        self.sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.logPosterior, backend=backend)
        self.sampler.run_mcmc(p0, self.iterations, progress=True)

        # get chain and log_prob in one-dimensional array (merged chains with burn-in)
        samples = self.sampler.get_chain(flat=True, discard=100)

        # get median mean and std
        self.posterior_mu = np.median(samples[:, 0])
        self.posterior_sigma = np.median(samples[:, 1])
        self.posterior_mu_std = np.std(samples[:, 0])
        return self.posterior_mu, self.posterior_sigma, self.posterior_mu_std

    ### Plotting figures ###

    def plot_emcee_chain(self, output_filename ='Chain_trace'):
        if not self.sampler:
            self.get_posterior_distribution()
        print('Plotting emcee chain trace to', output_filename)
        chain = self.sampler.get_chain()
        plt.close()
        fig, axs = plt.subplots(2)

        walker = 0
        chainT = chain.T
        while walker < self.nwalkers:
            axs[0].plot(chainT[0][walker], lw=0.5)
            axs[1].plot(chainT[1][walker], lw=0.5)
            walker += 1
        axs[0].set_ylabel('Mean')
        axs[1].set_ylabel('Std')
        fig.savefig(self.output_path + output_filename + '.pdf')
        fig.savefig(self.output_path + output_filename + '.png')
        plt.close()

    def plot_distribution(self, output_filename='Distribution_comparison'):
        print('Plotting distribution comparison to', output_filename)
        # Before BI
        fig = sns.kdeplot(self.y, label='Original distribution of HLs', color='#ff8c40', fill=True, alpha=.6,
                          linewidth=0)
        fig.axvspan(self.LOQ_lower, self.LOQ_upper, color='grey', alpha=0.3, lw=0)
        max_y = 0.8
        fig.plot([np.mean(self.y), np.mean(self.y)], [0, max_y], label='Original mean', color='#d3671b', ls='--')
        fig.plot([np.median(self.y), np.median(self.y)], [0, max_y], label='Original median', color='#d3671b',
                 ls=':')
        # After BI
        sns.kdeplot(np.random.normal(loc=self.posterior_mu, scale=self.posterior_sigma, size=500000),
                    label='Posterior distribution of HLs', color='#008e9b', fill=True, alpha=.5, linewidth=0)
        fig.plot([self.posterior_mu, self.posterior_mu], [0, max_y], label='Posterior mean',
                 color='#00727f', ls='--')
        # Plot actual data points
        y_list = np.random.uniform(0.01, 0.2, size=len(self.y))
        fig.scatter(x=self.y, y=y_list, color="black", marker='.')
        # print figure
        fig.set_xlabel(r"log(DT50)")
        fig.set_ylabel(r"p(log(DT50))")
        fig.legend()
        this = fig.get_figure()
        fig.savefig(self.output_path + output_filename + '.pdf')
        fig.savefig(self.output_path + output_filename + '.png')
        plt.close()

    def corner_plot(self, output_filename='Corner_plot'):
        import corner
        print("Plotting corner plot to {}".format(output_filename))
        assert self.save_backend == True, "ERROR: saved backend needed for this analysis: set_save_backend(True)"
        burnin = 100
        thin = 1
        samples = self.sampler.get_chain(discard=burnin, flat=True, thin=thin)
        log_prob_samples = self.sampler.get_log_prob(discard=burnin, flat=True, thin=thin)

        all_samples = np.concatenate((samples, log_prob_samples[:, None]), axis=1)

        labels = ["mean", "std", "log prob"]

        fig = corner.corner(all_samples, labels=labels)
        fig.savefig(self.output_path + output_filename + '.pdf')
        fig.savefig(self.output_path + output_filename + '.png')
        plt.close()

    def plot_prior_probabilities(self, output_filename=''):
        if output_filename == '':
            output_filename = 'Prior_probability_plot_{}_{}'.format(self.get_prior_mu(), self.get_prior_sigma())
        print("Plotting prior probability to {}".format(output_filename))
        range_theta = np.arange(-5, 7, 0.01)
        range_sigma = np.arange(0,2, 0.001)
        p_theta = []
        p_sigma = []
        [p_theta.append(np.exp(self.get_prior_probability_theta(theta))) for theta in range_theta]
        [p_sigma.append(np.exp(self.get_prior_probability_sigma(sigma))) for sigma in range_sigma]

        #create data frame
        df_theta = pd.DataFrame()
        df_theta[r'$\mu$'] = range_theta
        df_theta['Probability density'] = p_theta
        df_sigma = pd.DataFrame()
        df_sigma[r'$\sigma + \sigma_{min}$'] = range_sigma
        df_sigma['Probability density'] = p_sigma

        #plot
        fig, axes = plt.subplots(1, 2, figsize=(7, 3))
        sns.lineplot(data=df_theta, x=r'$\mu$', y='Probability density', color='#324b4e', ax=axes[0])
        sns.lineplot(data=df_sigma, x=r'$\sigma + \sigma_{min}$', y='Probability density', color='#008e9b', ax=axes[1])
        fig.tight_layout()
        fig.savefig(self.output_path + output_filename + '.pdf')
        fig.savefig(self.output_path + output_filename + '.png')
        plt.close()

    def plot_prior(self, output_filename=''):
        if output_filename == '':
            output_filename = 'Prior_plot_{}_{}'.format(self.get_prior_mu(), self.get_prior_sigma())
        print("Plotting prior distribution to {}".format(output_filename))

        fig, axes = plt.subplots(1, 2, figsize=(7, 3.3))
        mu_mean, mu_std = self.get_prior_mu()
        sigma_mean, sigma_std = self.get_prior_sigma()
        left = sns.kdeplot(np.random.normal(loc=mu_mean, scale=mu_std, size=50000),
                           label='{}={}, {}={}'.format(r'$\mu_{mean}$', mu_mean, r'$\mu_{std}$', mu_std), color='#324b4e', fill=True,
                           alpha=.8,
                           linewidth=0, ax=axes[0])
        left.set(title=r'$\mu$')
        left.axvspan(self.get_LOQ_lower(), self.get_LOQ_upper(), color='grey', alpha=0.3, lw=0,
                     label='LOQ range')
        sigma_dist = np.random.lognormal(mean=sigma_mean, sigma=sigma_std, size=50000)
        print('Stats: Median={}, Mean={}, Std={}'.format(np.median(sigma_dist), np.mean(sigma_dist),
                                                         np.std(sigma_dist)))
        sns.kdeplot(sigma_dist,
                            label='{}={}, {}={}'.format(r'$\sigma_{mean}$', sigma_mean,r'$\sigma_{std}$', sigma_std), color='#008e9b',
                            fill=True, alpha=.8, linewidth=0, ax=axes[1]).set(
            title=r'$\sigma + \sigma_{min}$')

        fig.legend(loc='upper left', bbox_to_anchor=(0.12, 0.89))
        fig.savefig(self.output_path + output_filename + '.pdf')
        fig.savefig(self.output_path + output_filename + '.png')
        plt.close()


### Utility functions ###
def get_normal_distribution(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


