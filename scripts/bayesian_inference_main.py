# This script was used to perform data analyses and create figures in Hafner et al., 2023
import sys
sys.path.insert(0, '../')
from pepper.pepper import Pepper
from pepper.datastructuresoil import DataStructureSoil
from pepper.bayesian import Bayesian

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import random

def figure_1(df_red, output_path):
    # mean, std
    color_prob = '#ff8c40'
    color_BI = '#008e9b'
    figure, axes = plt.subplot_mosaic([['left', 'right'], ['bottom', 'bottom']],
                                      constrained_layout=True)  # 1, 2, figsize=(7, 3))  # rows, columns
    sns.set_context("paper")
    sns.kdeplot(data=df_red, x='DT50_log_gmean', fill=True, common_norm=False, alpha=.7, linewidth=0,
                color=color_prob,
                cut=0, ax=axes['left'], label='Descriptive')
    sns.kdeplot(data=df_red, x='DT50_log_bayesian_mean', fill=True, common_norm=False, alpha=.7, linewidth=0,
                color=color_BI, cut=0, ax=axes['left'], label='Inferred').set(xlabel='Mean')
    sns.kdeplot(data=df_red, x='DT50_log_std', fill=True, common_norm=False, alpha=.7, linewidth=0,
                color=color_prob, cut=0, ax=axes['right'], label='Descriptive')
    sns.kdeplot(data=df_red, x='DT50_log_bayesian_std', fill=True, common_norm=False, alpha=.7, linewidth=0,
                color=color_BI, cut=0, ax=axes['right'], label='Inferred').set(xlabel='Standard deviation')
    axes['right'].legend()
    # boxplot
    # sns.set_style('whitegrid')
    # figure, axes = plt.subplots(figsize=(7, 3))  # rows, columns
    # sns.set_context("paper")
    color_prob = '#ff944d'
    color_BI = '#2f929c'
    D = {'Standard deviation': [], 'Number of half-life values': [], 'Method': []}
    df_red.sort_values(by=['DT50_count'], inplace=True)
    for index, row in df_red.iterrows():
        D['Number of half-life values'].append(get_category(row['DT50_count']))
        D['Number of half-life values'].append(get_category(row['DT50_count']))
        D['Method'].append('Descriptive')
        D['Method'].append('Inferred')
        if row['DT50_count'] < 3:
            D['Standard deviation'].append(np.NaN)
        else:
            D['Standard deviation'].append(row['DT50_log_std'])
        D['Standard deviation'].append(row['DT50_log_bayesian_std'])
    df = pd.DataFrame(D)
    sns.boxplot(data=df, x='Number of half-life values', y='Standard deviation', hue='Method',
                hue_order=['Descriptive', 'Inferred'],
                palette=sns.color_palette([color_prob, color_BI, '#ffffff'], 3, ),
                ax=axes['bottom'])
    # plt.legend()
    # plt.tight_layout()
    plt.savefig(output_path + '.pdf')
    plt.savefig(output_path + '.png')
    plt.close()

def figure_2(prior, output_path):
    analyse_prior(prior, output_path)

def figure_3(prior, output_path):
    # define compounds to visualize:
    # index	compound_id	smiles
    # 600	https://envipath.org/package/5882df9c-dae1-4d80-a40e-db4724271456/compound/48fcafd1-7b0c-4149-8d35-3f4cc1aab681	CN1CN(C)C(=S)SC1
    y1 = [-0.552841969, -0.26760624, 0.113943352, -3.537602002, -2.080921908, -1.22184875, -0.22184875, -3.537602002,
          -0.522878745]
    n1 = 'Dazomet'
    # index	compound_id	smiles
    # 845	https://envipath.org/package/5882df9c-dae1-4d80-a40e-db4724271456/compound/16a66ffb-d8e6-43fe-92c9-0d5a25b94821	CC1=NC(=NC(=N1)N)OC
    y2 = [2.04060234, 2.033825694, 2.06228107, 1.635483747, 2.28057837, 2.009450896, 1.804139432, 2.357934847,
          2.479863113, 1.718501689, 2.469822016, 1.81756537, 1.462397998, 2.491921713, 2.346607217, 1.81756537,
          2.522704993, 1.833784375, 1.718501689, 1.828015064, 1.804820679, 2.306425028, 2.564748892, 3, 1.919078092,
          2.055378331, 1.930439595, 2.155609284, 1.63748973, 1.791690649, 1.63748973, 1.462397998, 3, 2.005609445,
          2.515979756, 1.342422681, 2.396896449, 2.28057837, 3, 1.595496222, 3, 0.857332496, 1.342422681, 2.072984745,
          2.304705898, 2.009450896, 2.275080898, 2.009450896, 2.511482289, 2.415140352, 2.396896449, 2.07809415,
          2.072984745, 2.396896449, 2.005609445, 1.350248018, 1.804139432, 1.591064607, 1.650307523]
    n2 = 'IN-A4098 (Triazine amine)'
    # 1404 CC(C(=O)OCC1CCCO1)OC2=CC=C(C=C2)OC3=NC4=CC=C(C=C4N=C3)Cl
    y3 = np.log10([1, 1, 7.87, 1, 1, 1, 1, 1, 1, 1, 1])
    n3 = 'Quizalofop-P-tefuryl'
    # 2370	https://envipath.org/package/5882df9c-dae1-4d80-a40e-db4724271456/compound/e14efb46-3e84-453f-a96f-a0090b6441be	COC1=CC=NC(=C1O)C(=O)O
    y4 = [3, 3, 3, 3]
    n4 = 'X696476'
    # 2320	https://envipath.org/package/5882df9c-dae1-4d80-a40e-db4724271456/compound/1c15c524-6ab3-4889-81e9-33c8910d798a	CC1=NC(=NC(=N1)N(C)C(=O)NS(=O)(=O)C2=C(C=CC=C2)C(=O)OC)OC
    y5 = np.log10([14, 3, 5, 46, 5, 30, 4, 10, 12, 14, 5, 5, 10])
    n5 = 'Tribenuron-methyl'
    # 2608	https://envipath.org/package/5882df9c-dae1-4d80-a40e-db4724271456/compound/1ccac3c2-e8d1-4fa2-b11e-60011a8381b0	C[C@H]1CN(C[C@H](C)O1)C2CCCCCCCCCCC2
    y6 = np.log10([13, 10])  # only two values
    n6 = 'trans Dodemorph'
    # 2582	https://envipath.org/package/5882df9c-dae1-4d80-a40e-db4724271456/compound/5855469d-5893-4ff2-b804-76d589d7bc73	CCC(C)S(=O)(=O)C
    y7 = np.log10([4.5])  # only one value
    n7 = 'methyl-2-butyl sulfone'
    # 355 https://envipath.org/package/5882df9c-dae1-4d80-a40e-db4724271456/compound/000160cf-d55b-4b96-87b8-4c1b1193d425
    y8 = np.log10([1949])  # only one value
    n8 = 'Butralin'
    # 734	https://envipath.org/package/5882df9c-dae1-4d80-a40e-db4724271456/compound/64d2fb01-2cd9-415e-b604-6a416e50cc2e	CCCCOC(=O)[C@@H](C)OC1=CC=C(C=C1)OC2=CC=C(C=N2)C(F)(F)F
    y9 = np.log10([0.083, 0.083])
    n9 = 'Fluazifop-P-butyl '

    y_list = [y7, y6, y5, y1, y3, y2, y8, y9, y4]
    y_names = [n7, n6, n5, n1, n3, n2, n8, n9, n4]

    plot_multiple_distributions(y_list, y_names, prior, output_path=output_path)

def run_bayes_simple(y, prior, comment=[]):
    print_y_stats(y)
    bayesian = Bayesian(y=y, comment_list=comment)
    bayesian.set_save_backend(True)
    bayesian.set_iterations(2000)
    bayesian.set_prior_mu(mean=prior['mu_mean'], std=prior[
        'mu_std'])  # Can I categorize and use -1 if lower_LOQ = -1, 3 if upper_LOQ = 3, and 1.5 if no LOQ applies?
    bayesian.set_prior_sigma(mean=prior['sigma_mean'],
                             std=prior['sigma_std'])  # Shall I use different std priors for different data quality?
    bayesian.plot_emcee_chain()
    bayesian.plot_distribution()
    bayesian.plot_prior()
    # bayesian.log_prob_distribution_plot()
    bayesian.corner_plot()
    bayesian.plot_prior_probabilities()

    new_mean, new_std = bayesian.get_posterior_distribution()
    print('\nPrior mean: {}'.format(bayesian.get_prior_mu()))
    print('Prior std: {}'.format(bayesian.get_prior_sigma()))
    print('\nBayesian inference mean:', new_mean)
    print('Bayesian inference std:', new_std)

def plot_multiple_distributions(y_list, y_names, prior, y_comments=None, output_path='Case_studies_BI',
                                number_of_columns=3):
    size = len(y_list)
    if size % number_of_columns == 0:
        number_of_rows = size // number_of_columns
    else:
        number_of_rows = size // number_of_columns + 1
    figure, axes = plt.subplots(number_of_rows, number_of_columns,
                                figsize=(6 * number_of_columns, 2.5 * number_of_rows),
                                sharex=True)  # rows, columns #  figsize=(20, 20)
    row_counter = 0
    column_counter = 0
    for index, y in enumerate(y_list):
        compound_name = y_names[index]
        if y_comments:
            comment = y_comments[index]
        else:
            comment = []
        print_y_stats(y)
        ax = axes[row_counter, column_counter]
        bayesian = Bayesian(y=y, comment_list=comment)
        bayesian.set_iterations(2000)
        bayesian.set_prior_mu(mean=prior['mu_mean'], std=prior['mu_std'])
        bayesian.set_prior_sigma(mean=prior['sigma_mean'], std=prior['sigma_std'])
        bayesian.set_lower_limit_sigma(0.2)

        post_mean, post_std, post_mean_std = bayesian.get_posterior_distribution()
        print('\nBayesian inference mean:', post_mean)
        print('Bayesian inference std:', post_std)
        print('Bayesian inference mean_std:', post_mean_std)

        # get chain
        chain = bayesian.get_sampler().get_chain()
        chainT = chain.T
        walker = 0
        mean_list = []
        std_list = []
        while walker < bayesian.nwalkers:
            mean_list.extend(chainT[0][walker][bayesian.get_burnin():])
            std_list.extend(chainT[1][walker][bayesian.get_burnin():])
            walker += 1
        index_list = range(0, len(mean_list), 1)
        sample_indices = random.sample(index_list, 1000)

        # make figure
        title = "{}, n={}".format(compound_name, len(y))
        # fig = sns.kdeplot(y, label='Descriptive distribution', color='#ff8c40', fill=True, alpha=.6, linewidth=0, ax=ax)

        # After BI
        fig = sns.kdeplot(np.random.normal(loc=mean_list[0], scale=std_list[0], size=2000),
                          label='Sampled distributions', color='#2f929c', fill=False, alpha=.1, linewidth=0.01,
                          ax=ax)
        fig.set(title=title)
        for i in sample_indices[1:]:
            sns.kdeplot(np.random.normal(loc=mean_list[i], scale=std_list[i], size=2000),
                        color='#008e9b', fill=False, alpha=.3, linewidth=0.1, ax=ax)

        fig.axvspan(bayesian.get_LOQ_lower(), bayesian.get_LOQ_upper(), color='grey', alpha=0.3, lw=0,
                    label='LOQ range')

        _, max_y = fig.get_ylim()
        min_x, _ = fig.get_xlim()

        fig.plot([np.mean(y), np.mean(y)], [0, max_y], label='Descriptive mean', color='#d3671b', ls='--')
        fig.plot([np.median(y), np.median(y)], [0, max_y], label='Descriptive median', color='#d3671b', ls=':')

        fig.plot([post_mean, post_mean], [0, max_y], label='Inferred mean',
                 color='#00727f', ls='--')
        # Plot actual data points
        real_y = []
        beyond_y = []
        for val in y:
            if val <= bayesian.get_LOQ_lower():
                beyond_y.append(val)
            elif val >= bayesian.get_LOQ_upper():
                beyond_y.append(val)
            else:
                real_y.append(val)
        y_list_real = np.random.uniform(0.01, 0.2, size=len(real_y))
        y_list_beyond = np.random.uniform(0.01, 0.2, size=len(beyond_y))
        # plot data points within LOQ range
        fig.plot(real_y, y_list_real, color="black", marker='.', linewidth=0, label='Valid data points')
        # plot data points beyond LOQ
        fig.plot(beyond_y, y_list_beyond, color="black", marker='x', linewidth=0, label='Censored data points')
        # print figure
        fig.set_xlabel(r"log(DT50)")
        fig.set_ylabel(r"p(log(DT50))")
        plt.xlim([-5, 6])
        if len(y) >= 3:
            original_std = np.round(np.std(y), 2)
        else:
            original_std = 'ND'
        fig.annotate('$\mathbf{{Descriptive:}}$\nmean={}\nstd={}'.format(np.round(np.mean(y), 2), original_std) +
                     '\n$\mathbf{{Inferred:}}$\n{}={}\n{}={}\n{}={}'.format(r'$\mu_{mean}$', np.round(post_mean, 2),
                                                                            r'$\mu_{std}$',
                                                                            np.round(post_mean_std, 2),
                                                                            r'$\sigma_{mean}$',
                                                                            np.round(post_std, 2)),
                     xy=(-4.8, max_y - (0.65 * max_y)))
        if column_counter == number_of_columns - 1:
            column_counter = 0
            row_counter += 1
        else:
            column_counter += 1

    plt.legend(loc='lower right', bbox_to_anchor=(0.88, 0), bbox_transform=figure.transFigure, ncol=8)
    plt.savefig(output_path + '.pdf')
    plt.savefig(output_path + '.png', dpi = 600)

def print_y_stats(y):
    print('Descriptive mean:', np.mean(y))
    print('Descriptive median:', np.median(y))
    print('Descriptive std:', np.std(y))

def analyse_prior(prior, output_path):
    y = [0, 1]
    comment = ['', '']
    bayesian = Bayesian(y=y, comment_list=comment)
    bayesian.set_prior_mu(mean=prior['mu_mean'], std=prior[
        'mu_std'])  # Can I categorize and use -1 if lower_LOQ = -1, 3 if upper_LOQ = 3, and 1.5 if no LOQ applies?
    bayesian.set_prior_sigma(mean=prior['sigma_mean'],
                             std=prior['sigma_std'])  # Shall I use different std priors for different data quality?
    bayesian.set_lower_limit_sigma(0.2)
    bayesian.plot_prior_probabilities(output_filename = output_path)

def get_category(hl_count):
    if 1 <= hl_count <= 8:
        category = str(hl_count)
    elif 9 <= hl_count <= 10 :
        category = '9-10'
    elif 11 <= hl_count <= 15 :
        category = '11-15'
    elif 16 <= hl_count <= 20:
        category = '16-20'
    elif 21 <= hl_count <= 60:
        category = '21-60'
    else:
        raise ValueError('No matching category for {}'.format(hl_count))
    return category

### MAIN ###
if __name__ == '__main__':
    # initiate
    pep = Pepper()
    pep.set_tag('all_data')

    # set important column names
    pep.set_target_variable_name('logDT50_mean')
    pep.set_target_variable_std_name('logDT50_std')
    pep.set_smiles_name('SMILES')
    pep.set_id_name('ID')

    # load data
    soil_data = DataStructureSoil(pep)
    soil_data.curate_annotate(from_csv=True, from_paper=True)
    soil_data.reduce_data()

    # Reproduce figures from paper
    pep.set_data_directory('/pepper_data/bayesian_inference/')
    figure_1(soil_data.cpd_data, pep.get_data_directory() + 'figure_1')

    # define prior
    prior = {
        'mu_mean': 1,
        'mu_std': 2,
        'sigma_mean': 0.4,
        'sigma_std': 0.4
    }

    figure_2(prior, pep.get_data_directory() + 'figure_2')
    figure_3(prior, pep.get_data_directory() + 'figure_3')