import os
import numpy as np
from copy import deepcopy
import pandas as pd
import matplotlib
import pickle

from matplotlib.pyplot import tight_layout
from scipy.ndimage import label

from pepper_lab.pepper import Pepper
from sklearn.metrics import mean_squared_error, r2_score
from pepper_lab.util import *
# from openTSNE.tsne import TSNE

# # Umap related imports
# import umap
# import umap.plot
from sklearn.cluster import KMeans
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS

matplotlib.use('Agg')


class Visualize(Pepper):
    def __init__(self, visualized_object, analysis_type):
        """Section could be data object, descriptors, model, modelling"""
        super().__init__()

        self.object = visualized_object
        self.object_name = visualized_object.get_object_name()  # the data structure to be visualized
        self.set_data_directory(os.path.join(visualized_object.pepper.data_directory,
                                             self.object_name, 'visualization', analysis_type))

        # used to specify output file name
        self.tag = visualized_object.get_tag()  # Changed 'object.pepper' to 'object' to control the tag based on modeling
        self.data_type = visualized_object.get_data_type()
        self.setup_name = visualized_object.get_setup_name()

        #output settings
        self.context = 'paper'  # other option: presentation
        self.color_palette_3 = ["#428aa3", "#c3618b", "#7a8838"]
        self.color_palette_gradient = ["#428aa3", "#4b87b0", "#6581b8", "#8579b6", "#a46eaa", "#b96696",
                                       "#c5627e", "#c86365", "#bf6b4c", "#ae753a", "#968032", "#7a8838"]
        self.color_palette_models = ['#5ab7e1', '#0082aa', '#fff189', '#ae7a0d', '#9fe1d5', '#003c34', '#f380f0', '#5e0062', '#d4a375', '#9b6f44', '#d5cabd', '#50434f']
        self.color_palette_descriptors = ['#f380f0', '#810082', '#74b8ac', '#005249', '#a2acbd', '#2c4a6e', '#ffedcb', '#dda11d']
        self.color_palette_data_sets = ['#60baae', '#2f4858', '#006b61', '#2f4858', '#92a19f', '#554516', '#b49205']

        # openTSNE related
        self.embedding_test = pd.DataFrame()
        self.embedding_for_plot = pd.DataFrame()

    @staticmethod
    def fetch_object(self):
        pass

    def set_context(self):
        if type in ['paper', 'presentation']:
            self.context = type
        else:
            return ValueError

    def set_setup_name(self, setup_name: str):
        self.setup_name = setup_name

    def modeling_summary(self):
        training_directory = self.object_path  # todo:
        df = pd.read_csv(training_directory, sep='\t')


    def load_df_from_path(self, file_name: str): #todo: what is this function for?
        object_df = pd.read_csv(self.object_path + file_name)

    def scatterplot_predicted_vs_test(self):
        assert self.object_name == 'Model', \
            "This function cannot be applied to the object {}".format(self.object_name)
        sns.set(rc={"figure.figsize": (5, 5)})
        sns.set_theme(style="whitegrid")
        sns.set_style("ticks")
        sns.set_context(self.context)
        output_file_path = os.path.join(self.get_data_directory(),
                                 'scatterplot_pred_vs_test_{}_{}_{}_{}.pdf'.format(self.data_type, self.tag,
                                                                              self.setup_name, self.object.settings_string))

        assert len(self.object.y_pred) != 0, 'No prediction available to be visualized for this model'
        assert len(self.object.y_test) == len(self.object.y_pred), 'The predicted and the true y have not the same size'

        df = pd.DataFrame()
        df['Experimental'] = self.object.y_test
        df['Predicted'] = self.object.y_pred
        if len(self.object.y_pred_score) != 0 :
            print('using prediction score')
            df['Prediction_score'] = self.object.y_pred_score
            ax = sns.scatterplot(x='Experimental', y='Predicted', hue='Prediction_score', data=df, palette="coolwarm")
            ax.errorbar(x=df['Experimental'],y=df['Predicted'],yerr=df['Prediction_score'], fmt='.',ms=0, alpha=0.1, elinewidth=0.5, color='black')
        else:
            ax = sns.scatterplot(x='Experimental', y='Predicted', data=df)
        if self.object.target_variable_std_name:
            df['Experimental_std'] = self.object.y_test_std
            df['Experimental_std'] = df['Experimental_std'].fillna(0)
            try:
                print('trying to add error bars')
                ax.errorbar(x=df['Experimental'], y=df['Predicted'], xerr=df['Experimental_std'], fmt='.', ms=0, alpha=0.1, elinewidth=0.5, color='black')
            except ValueError:
                print('Are there nan values in the the standard deviation?')
                print(self.object.y_test_std.isna().sum())
                pass

        # Add +/-1 log unit lines
        # min_val = min(df.min().values)
        # min_val = df.Experimental.min()
        min_val = np.array(df).min()*1.1
        # max_val = max(df.max().values)
        # max_val = df.Experimental.max()
        max_val = np.array(df).max()*1.1
        ax.plot([min_val, max_val], [min_val, max_val], ls='-', color='black')
        ax.plot([min_val + 1, max_val], [min_val, max_val - 1], ls='--', color='black')
        ax.plot([min_val, max_val - 1], [min_val + 1, max_val], ls='--', color='black')
        ax.set(xlim=[min_val, max_val], ylim=[min_val, max_val])

        print('Saving figure to', output_file_path)
        plt.savefig(output_file_path)
        plt.close()

    def boxplot_scatterplot_CV_performance(self, by_categories=['regressor', 'descriptors']):
        """
        Plot cross-validation performance by different categories in a scatter_box_plot
        @param by_categories: 'regressor', 'descriptors', 'train_fraction'
        """
        assert self.object_name == 'Modeling', \
            "This function cannot be applied to the object {}".format(self.object_name)
        data = self.object.test_scores
        data.rename(columns={'r2': 'R2', 'rmse': 'RMSE'}, inplace =True)
        for category in by_categories:
            self.scatter_box_plot(data, category, self.get_color_palette_by_category(category))


    def boxplot_CV_performance(self, categories: list):
        """
        Plot cross-validation performance for regressors x descriptors
        @param categories: 'regressor', 'descriptors', 'train_fraction'
        """
        assert self.object_name == 'Modeling', \
            "This function cannot be applied to the object {}".format(self.object_name)
        data = self.object.test_scores
        data.rename(columns={'r2': 'R2', 'rmse': 'RMSE'}, inplace =True)

        self.performance_box_plots(data=data, categories=categories, palette=self.get_color_palette_by_category(categories[0]))

    def get_color_palette_by_category(self,category:str):
        if category == 'regressor':
            pal = self.color_palette_models
        elif category == 'descriptors':
            pal = self.color_palette_descriptors
        elif category == 'feature_selection':
            pal = self.color_palette_gradient
        elif category == 'train_fraction':
            pal = self.color_palette_gradient
        else:
            pal = self.color_palette_gradient
        return pal

    def plot_experimental_performance_simulation(self, df):
        figure, axes = plt.subplots(1, 2, figsize=(5, 4))  # rows, columns
        sns.set_context('paper')
        df_R2 = df.loc[:, ['R2_exp', 'R2_dist']]
        df_R2.rename(columns={'R2_exp': 'from\nreported\nvalues', 'R2_dist': 'from\ndistribution'}, inplace=True)
        ax = sns.boxplot(data=df_R2, ax=axes[0], palette='Paired').set_title('R2')
        df_RMSE = df.loc[:, ['RMSE_exp', 'RMSE_dist']]
        df_RMSE.rename(columns={'RMSE_exp': 'from\nreported\nvalues', 'RMSE_dist': 'from\ndistribution'}, inplace=True)
        ax = sns.boxplot(data=df_RMSE, ax=axes[1], palette='Paired').set_title('RMSE')
        axes[0].set_ylim(0.6, 1)
        axes[1].set_ylim(0.0, 0.35)
        figure.tight_layout()
        output_filename = os.path.join(self.get_data_directory(),
                                 'estimation_possible_performance_{}_{}_{}.pdf'.format(self.data_type, self.tag,
                                                                              self.setup_name))
        print("Saving figure to {}".format(output_filename))
        plt.savefig(output_filename)
        plt.close()

    # distribution of target variable
    def plot_target_variable_distribution(self, mean_name: str, std_name: str,
                                          cutoff_value: int = 10,
                                          include_BI: bool = False,
                                          BI_mean_name: str = 'DT50_log_bayesian_mean',
                                          BI_std_name: str = 'DT50_log_bayesian_std'):
        """
        This function plots the distribution of the target variable(s), for the whole data set and for a reduced data set of compounds withmore than k (cutoff_value) reported target variables.
        @param mean_name: The column name where mean (gmean, median, etc) values are reported. The mean is obtained for a same compound, but different experiments.
        @param std_name: The column name where the standard deviation of the target variable is reported
        @param cutoff_value: The minimum number of data points to be used for the calculation of a realistic standard deviation
        @param include_BI: Whether Bayesian inference has been performed already or not - if yes, also provide BI_mean_name and BI_std_name
        @param BI_mean_name: The column where the Bayesian-inferred mean can be found
        @param BI_std_name: The column where the Bayesian-inferred standard deviation can be found
        """
        figure, axes = plt.subplot_mosaic([['left', 'right']],
                                          constrained_layout=True, figsize=(7, 3))  # 1, 2, figsize=(7, 3))  # rows, columns
        sns.set_context("paper")
        df = self.object.cpd_data
        df_top_k = self.object.cpd_data.loc[self.object.cpd_data["DT50_count"]>=cutoff_value]
        print("Dataset with {} or more reported experimental values:".format(cutoff_value))
        print("\t - Number of compounds: {}".format(len(df_top_k)))
        print("\t - Distribution of the target variable: mean = {}, std = {}".format(
            np.round(np.mean(df_top_k[mean_name]), 2),
            np.round(np.std(df_top_k[mean_name]), 2)))
        print("\t - Distribution of the standard deviation of the target variable: mean = {}, std = {}".format(
            np.round(np.mean(df_top_k[std_name]), 2),
            np.round(np.std(df_top_k[std_name]), 2)))
        sns.kdeplot(data=df, x=mean_name, fill=True, common_norm=False, alpha=.6, linewidth=0,
                    color=self.color_palette_3[0],
                    cut=0, ax=axes['left'], label='Descriptive')
        sns.kdeplot(data=df_top_k, x=mean_name, fill=True, common_norm=False, alpha=.6, linewidth=0,
                    color=self.color_palette_3[2], cut=0, ax=axes['left'], label='Descriptive (n>={})'.format(cutoff_value)).set(xlabel='Mean')
        sns.kdeplot(data=df, x=std_name, fill=True, common_norm=False, alpha=.6, linewidth=0,
                    color=self.color_palette_3[0], cut=0, ax=axes['right'], label='Descriptive')
        sns.kdeplot(data=df_top_k, x=std_name, fill=True, common_norm=False, alpha=.6, linewidth=0,
                    color=self.color_palette_3[2], cut=0, ax=axes['right'], label='Descriptive (n>={})'.format(cutoff_value)).set(xlabel='Standard deviation')
        if include_BI:
            sns.kdeplot(data=df, x=BI_mean_name, fill=True, common_norm=False, alpha=.6, linewidth=0,
                    color=self.color_palette_3[1], cut=0, ax=axes['left'], label='Inferred')
            sns.kdeplot(data=df, x=BI_std_name, fill=True, common_norm=False, alpha=.6, linewidth=0,
                        color=self.color_palette_3[1], cut=0, ax=axes['right'], label='Inferred')
        axes['right'].legend()
        output_filename = os.path.join(self.get_data_directory(),
                                 'plot_target_variable_distribution_{}_{}_{}.pdf'.format(self.data_type, self.tag,
                                                                              self.setup_name))
        print("Saving figure to {}".format(output_filename))
        plt.savefig(output_filename)
        plt.close()

    def scatter_box_plot(self,data: pd.DataFrame, category: str ,palette: list):
        """
        Visualize performance with a scatter plot RMSE vs R2, R2 box plot, and RMSE box plot colored by category
        @param data: data frame with performance scores
        @param category: category to be visualized: 'regressor', 'descriptors', or 'train_fraction'
        @param palette: color palette
        """
        sns.set(rc={"figure.figsize": (5, 5)})
        sns.set_theme(style="whitegrid")
        sns.set_style("ticks")
        sns.set_context(self.context)

        # adjust color palettes
        n_colors_needed = data.groupby(category).size().size
        assert len(palette)>= n_colors_needed, "Palette has not enough colors: {}".format(palette)
        palette = palette[:n_colors_needed]

        # plot properties
        boxprops = {'edgecolor': 'k', 'linewidth': 1}
        lineprops = {'color': 'k', 'linewidth': 1}

        boxplot_kwargs = {'boxprops': boxprops, 'medianprops': lineprops,
                          'whiskerprops': lineprops, 'capprops': lineprops,
                          'width': 0.75, 'palette': palette, 'saturation':0.5}

        stripplot_kwargs = {'linewidth': 0.6, 'size': 6, 'alpha': 0.7,
                            'palette': palette}

        # plot figure
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        sns.scatterplot(data=data, y='RMSE', x='R2', hue=category, ax=axs[0], palette=palette)
        ax1 = sns.boxplot(y='RMSE', x=category, data=data, hue=category, ax=axs[1], **boxplot_kwargs)
        ax1 = sns.stripplot(y='RMSE', x=category, data=data, hue=category, ax=axs[1], **stripplot_kwargs)
        ax2 = sns.boxplot(y='R2', x=category, data=data, hue=category, ax=axs[2], **boxplot_kwargs)
        ax2 = sns.stripplot(y='R2', x=category, data=data, hue=category, ax=axs[2], **stripplot_kwargs)
        [ax.legend_.remove() for ax in [ax1, ax2]]
        # [ax.set_xticklabels(ax.get_xticklabels(), rotation=40) for ax in [ax1, ax2]]
        plt.tight_layout()
        # save file
        output_filename = os.path.join(self.get_data_directory(),
                                 'scatter_box_plot_{}_{}_{}_{}.pdf'.format(category, self.data_type, self.tag,
                                                                                       self.setup_name))
        print("Save figure to {}".format(output_filename))
        plt.savefig(output_filename)
        plt.close()

    def performance_box_plots(self, data: pd.DataFrame, categories: list, palette: list):
        """
        Visualize performance with a scatter plot RMSE vs R2, R2 box plot, and RMSE box plot colored by first category
        provided in parameter "categories"
        @param data: data frame with performance scores
        @param categories: list of category to be visualized. Possible items: 'regressor', 'descriptors',
                            'feature_selection', 'train_fraction'
        @param palette: color palette with enough colors to cover all values available for the first provided category (categories[0])
        """
        # create new column in data frame for combinations
        data['combination'] = data[categories].agg('\nx\n'.join, axis=1)
        number_of_boxes = len(data['combination'].unique())

        # adjust color palettes to order by categories[0]
        n_colors_needed = data.groupby(categories[0]).size().size
        assert len(palette)>= n_colors_needed, "Palette has not enough colors: {}".format(palette)
        palette = palette[:n_colors_needed]

        # replace names with display names
        replace_dict = {}
        for c in categories:
            names = data[c].unique()
            for this_name in names:
                replace_dict[this_name] = Util.get_display_name(this_name)

        # make names ready for display
        data.replace(replace_dict, inplace=True)
        # create new column in data frame for combinations
        data['combination'] = data[categories].agg('\nx\n'.join, axis=1)

        # set up plot
        sns.set_theme(style="whitegrid")
        sns.set_style("ticks")
        sns.set_context(self.context)


        boxprops = {'edgecolor': 'k', 'linewidth': 1}
        lineprops = {'color': 'k', 'linewidth': 1}

        boxplot_kwargs = {'boxprops': boxprops, 'medianprops': lineprops,
                          'whiskerprops': lineprops, 'capprops': lineprops,
                          'width': 0.75, 'palette': palette, 'saturation':0.5}

        stripplot_kwargs = {'linewidth': 0.6, 'size': 6, 'alpha': 0.7,
                            'palette': palette}

        # define figure size
        if number_of_boxes <16:
            figure_size = (8, 5)
        else:
            figure_size = (number_of_boxes*0.5, number_of_boxes*0.3)

        # plot
        fig, axs = plt.subplots(2, 1, figsize=figure_size)
        ax1 = sns.boxplot(y='RMSE', x='combination', data=data, hue = categories[0], ax=axs[0], **boxplot_kwargs)
        sns.stripplot(y='RMSE', x='combination', data=data, hue = categories[0], ax=axs[0], **stripplot_kwargs)
        ax2 = sns.boxplot(y='R2', x='combination', data=data, hue = categories[0], ax=axs[1], **boxplot_kwargs)
        sns.stripplot(y='R2', x='combination', data=data,hue = categories[0], ax=axs[1], **stripplot_kwargs)
        ax1.set_xticklabels([]) # do not show xtick labels on upper plot
        # ax2.set_xticklabels(ax2.get_xticklabels(), rotation=40) # can be used to rotate labels
        [ax.legend_.remove() for ax in [ax1, ax2]]
        [ax.axes.get_xaxis().set_label_text('') for ax in [ax1, ax2]]
        plt.tight_layout()
        output_filename = os.path.join(self.get_data_directory(),
                                 'performance_box_plots_{}_{}_{}_{}.pdf'.format('x'.join(categories), self.data_type, self.tag,
                                                                           self.setup_name))
        print("Save figure to {}".format(output_filename))
        plt.savefig(output_filename)
        plt.close()

    def feature_selection_plot(self):
        """
        Plots a boxplot of importance scores obtained from 5-fold CV for the selected features
        """
        assert self.object_name == 'Model', \
            "This function cannot be applied to the object {}".format(self.object_name)
        fig, axs = plt.subplots(figsize=(10, 5))
        sns.set_context('paper')
        # create color map by feature origin
        colors = {}
        data = self.object.selected_features_scores

        for index, row in data.iterrows():
            feature = row['Feature']
            colors[feature] = self.get_feature_space_color(self.object.descriptors.feature_space_map[feature])

        # settings
        # plot properties
        boxprops = {'edgecolor': 'k', 'linewidth': 1}
        lineprops = {'color': 'k', 'linewidth': 1}

        boxplot_kwargs = {'boxprops': boxprops, 'medianprops': lineprops, 'whiskerprops': lineprops,
                          'capprops': lineprops, 'width': 0.75, 'palette': colors, 'saturation':0.5}

        stripplot_kwargs = {'linewidth': 0.6, 'size': 6, 'alpha': 0.7, 'palette': colors}

        ax = sns.boxplot(data=data, x='Feature', y='Score', hue='Feature', **boxplot_kwargs)
        ax = sns.stripplot(data=data, x='Feature', y='Score', hue='Feature', **stripplot_kwargs)
        ax.set_xticks(ax.get_xticks())
        ax.set_xticklabels(ax.get_xticklabels(), rotation='vertical')

        # plot & save
        plt.tight_layout()
        output_filename = os.path.join(self.get_data_directory(),
                                 'feature_selection_plot_{}_{}_{}_{}.pdf'.format(self.data_type, self.tag,
                                                                           self.setup_name, self.object.settings_string))
        print("Save figure to {}".format(output_filename))
        plt.savefig(output_filename)
        plt.close()

    def plot_performance_vs_score_threshold(self, threshold_list = [0.5, 0.65, 0.8, 1], reverse = False):
        """
        Plot performance evolution as a function of on prediction score cutoff
        """
        assert self.object_name == 'Model', \
            "This function cannot be applied to the object {}".format(self.object_name)

        sns.set(rc={"figure.figsize": (5, 5)})
        sns.set_theme(style="whitegrid")
        sns.set_style("ticks")
        sns.set_context(self.context)
        df = self.object.predicted_target_variable
        assert df.get('predicted_score') is not None, "A prediction score must be provided"
        df.sort_values(by='predicted_score', inplace=True)
        colors = sns.color_palette('magma', 4)

        rmse_list = []
        r2_list = []
        cat_list = []
        std_list = []

        # iterate through list of predicted standard deviation thresholds to be evaluated
        for i in np.arange(start=0.0, stop=1, step=0.01):
            if reverse:
                this_df = df.loc[df['predicted_score'] >= i]
            else:
                this_df = df.loc[df['predicted_score'] <= i]
            y_pred = this_df['predicted']
            y_true = this_df['experimental']

            r2, rmse = self.get_scores(y_true, y_pred)

            rmse_list.append(rmse)
            r2_list.append(r2)
            cat_list.append(self.get_category(i, threshold_list, reverse))
            std_list.append(i)

        df_eval = pd.DataFrame()
        df_eval['R2'] = r2_list
        df_eval['RMSE'] = rmse_list
        df_eval['Threshold category'] = cat_list
        df_eval['Predicted std threshold'] = std_list

        fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))
        sns.set_theme("paper")
        sns.set_style("white")
        sns.scatterplot(data=df_eval, y='RMSE', x='Predicted std threshold', hue='Threshold category', palette=colors,
                        ax=axs[0])
        sns.scatterplot(data=df_eval, y='R2', x='Predicted std threshold', hue='Threshold category', palette=colors,
                        ax=axs[1])

        # plot & save
        plt.tight_layout()
        output_filename = os.path.join(self.get_data_directory(),
                                 'performance_vs_score_threshold_{}_{}_{}_{}.pdf'.format(self.data_type, self.tag,
                                                                              self.setup_name, self.object.settings_string))
        print("Save figure to {}".format(output_filename))
        plt.savefig(output_filename)
        plt.close()

    def scatterplots_by_thresholds(self, threshold_list = [0.5, 0.65, 0.8, 1], reverse = False):
        """
        Plot performance evolution depending on prediction score cutoff
        """
        assert self.object_name == 'Model', \
            "This function cannot be applied to the object {}".format(self.object_name)
        df = self.object.predicted_target_variable
        assert df.get('predicted_score') is not None, "A prediction score must be provided"
        output_filename = os.path.join(self.get_data_directory(),
                                 'scatterplots_by_threshold_{}_{}_{}_{}.pdf'.format(self.data_type, self.tag,
                                                                              self.setup_name, self.object.settings_string))
        fig, ax = plt.subplots(2, 2, figsize=(7, 7))
        sns.set_context("paper")
        sns.set_style("white")
        min_val = -3
        max_val = 4.5

        i = 0
        for t in threshold_list:
            this_df = deepcopy(df.loc[df['predicted_score'] <= t])

            this_df.rename(columns={'experimental': 'Experimental half-life [log(d)]', 'predicted': 'Predicted half-life [log(d)]'},
                           inplace=True)
            this_color = sns.color_palette('coolwarm', 13)[i*3+i]
            axis = sns.scatterplot(x='Experimental half-life [log(d)]', y='Predicted half-life [log(d)]',
                                   data=this_df, ax=self.get_axis(ax, i),
                                   color=this_color)
            axis.plot([min_val, max_val], [min_val, max_val], ls='-', color='black')
            axis.plot([min_val + 1, max_val], [min_val, max_val - 1], ls='--', color='black')
            axis.plot([min_val, max_val - 1], [min_val + 1, max_val], ls='--', color='black')


            y_true = this_df['Experimental half-life [log(d)]']
            y_pred = this_df['Predicted half-life [log(d)]']
            r2, rmse = self.get_scores(y_true, y_pred)
            axis.text(-3, 4.2, 'Threshold: {}'.format(t))
            axis.text(-3, 3.7, '# Compounds: {} ({}%)'.format(len(y_pred), round(len(y_pred)/len(df['predicted_score'])*100)))
            axis.text(-3, 3.2, 'R2: {}'.format(round(r2, 2)))
            axis.text(-3, 2.7, 'RMSE: {}'.format(round(rmse, 2)))

            axis.errorbar(x=this_df['Experimental half-life [log(d)]'], y=this_df['Predicted half-life [log(d)]'],
                         xerr=this_df['experimental_std'], yerr=this_df['predicted_score'], fmt='.', alpha=0.8,
                          color=this_color, elinewidth = 0.5)
            plt.subplots_adjust(hspace=0.3)
            i += 1
        # plt.tight_layout()
        print("Save figure to {}".format(output_filename))
        plt.savefig(output_filename)
        plt.close()

    @staticmethod
    def get_category(score, thresholds, reverse = False):
        for t in thresholds:
            if not reverse and score <= t:
                return f"<{t}"
            if reverse and score >= t:
                return f">{t}"

    @staticmethod
    def get_scores(y_true, y_pred):
        if len(y_pred) < 3 :
            return np.nan, np.nan
        else:
            r2 = r2_score(y_true, y_pred)
            rmse = mean_squared_error(y_true, y_pred)
            return r2, rmse

    @staticmethod
    def get_axis(ax, i):
        if i == 0:
            this_ax = ax[0][0]
        if i == 1:
            this_ax = ax[0][1]
        if i == 2:
            this_ax = ax[1][0]
        if i == 3:
            this_ax = ax[1][1]
        return this_ax

    @staticmethod
    def get_feature_space_color(feature):  #todo: @Jasmin: is this supposed to be static
        """
        Get a predefined color for each feature space
        @param feature: name of descriptor
        @return:
        """
        if feature == 'maccs':
            return '#aba300'
        elif feature == 'padel':
            return '#005b48'
        elif feature == 'ep_trig':
            return '#ff6f91'
        elif feature == 'ep_prob':
            return '#746e00'
        elif feature == 'mfps':
            return '#bea5a9'
        elif feature == 'mordred':
            return '#a00039'
        elif feature == 'clogp':
            return '#f7f5dd'
        elif feature == 'plant_fp':
            return '#ff655d'
        elif feature == 'rdkitfps':
            return '#32222b'
        elif feature[:2] == 'PC':
            return '#32222b'
        else:
            raise NotImplementedError(f'No color defined for {feature} feature space')

    # #------------------------------------#
    # # openTSNE related stuff  #
    # #------------------------------------#
    # def train_my_openTSNE(self, training_fingerprint=None, training_fingerprint_directory='', load_from_csv=False):
    #     if load_from_csv:
    #         training_fingerprint = pd.read_csv(training_fingerprint_directory)
    #
    #     tsne = TSNE(
    #         perplexity=100,
    #         n_iter=2000,
    #         metric='jaccard',
    #         random_state=42,
    #         verbose=True,
    #     )
    #     self.embedding_train = tsne.fit(training_fingerprint)
    #
    # def get_openTSNE_embedding(self, my_mfps, load_from_csv=False, load_embedding_from='directory.sav'):
    #     if load_from_csv:
    #         # load the model from disk
    #         filename = 'open_tsne_trained.sav'
    #         file_directory = os.path.join(load_embedding_from, filename)
    #         loaded_embedding = pickle.load(open(file_directory, 'rb'))
    #         self.embedding_train = loaded_embedding
    #
    #     # Get the morgan fingerprints but remove the smiles if present
    #     if self.object.smiles_name in my_mfps.columns:
    #         my_mfps.drop(self.object.smiles_name, axis=1, inplace=True)
    #
    #     # Transform the mfps
    #     self.embedding_test = self.embedding_train.transform(my_mfps)
    #
    #     self.embedding_for_plot['tsne_v1'] = self.embedding_test[:, 0]
    #     self.embedding_for_plot['tsne_v2'] = self.embedding_test[:, 1]

    def get_embedding_plot(self, plot_name='opentsne_embedding'):
        figure_directory = self.get_data_directory()
        ax = sns.scatterplot(data=self.embedding_for_plot, x='tsne_v1', y='tsne_v2', s=1, hue='class',
                             alpha=0.6, edgecolor='black')
        ax.legend(loc='upper left', bbox_to_anchor=(1.00, 0.75), ncol=1)
        plt.savefig(str(figure_directory) + '/{}.pdf'.format(plot_name), bbox_inches='tight', dpi=1200)

    # def show_chemical_space(self, my_mfps, plot_name, load_from_csv=False,  load_embedding_from='embedding_path'):
    #     self.get_openTSNE_embedding(my_mfps, load_from_csv)
    #     self.get_embedding_plot(plot_name=plot_name)

    # #------------------------------------#
    # # UMAP related stuff  #
    # #------------------------------------#
    #
    # def get_umap_plot(self):
    #     my_umap, umap_embedding = self.get_embedding()
    #     ax = sns.scatterplot(x=umap_embedding[:, 0], y=umap_embedding[:, 1])
    #     ax.legend(loc='upper left', bbox_to_anchor=(1.00, 0.75), ncol=1)
    #     plt.show()
    #
    # # Get the UMAP embedding
    # def get_embedding(self):
    #     my_umap = umap.UMAP(random_state=42,
    #                         n_neighbors=50,
    #                         min_dist=0.01,
    #                         n_components=2)
    #     umap_embedding = my_umap.fit_transform(self.object.features)
    #     return my_umap, umap_embedding
    #
    # # Display the diagnostic embedding
    #
    # def show_diagnose(self, my_umap):
    #     mapper = my_umap.fit(self.object.features)
    #     umap.plot.diagnostic(mapper, diagnostic_type='pca')
    #     plt.show()
    #
    # # Use the UMAP embedding is input in KNN(n=10) to get clusters
    # def get_knn_clusters(self, umap_embedding):
    #     kmeans_labels = KMeans(n_clusters=10).fit_predict(umap_embedding)
    #     my_kmeans_labels = ['C_' + str(x) for x in kmeans_labels]
    #     self.object.data['kmeans_labels'] = my_kmeans_labels
    #     return self.object.data, my_kmeans_labels
    #
    # # plot UMAP embedding with the clusters as colors
    # @staticmethod
    # def show_clusters(umap_embedding, my_kmeans_labels):
    #     ax = sns.scatterplot(x=umap_embedding[:, 0], y=umap_embedding[:, 1], hue=my_kmeans_labels, palette='viridis')
    #     ax.legend(loc='upper left', bbox_to_anchor=(1.00, 0.75), ncol=1)
    #     plt.show()

    # Print number of molecules in each cluster
    def get_clusters_dict(self):
        clusters_dict = {}
        for cluster_label in list(self.object.data.kmeans_labels.unique()):
            cluster_label_df = self.object.data[self.object.data.kmeans_labels == cluster_label]
            print("{}:{}".format(cluster_label, len(cluster_label_df.CanonicalSMILES.unique())))
            clusters_dict[cluster_label] = cluster_label_df
        return clusters_dict

    # For each cluster get and display the MCS (Maximum Common Substructure)
    @staticmethod
    def display_MCS(clusters_dict):
        for cluster_df in clusters_dict.items():
            mols = [Chem.MolFromSmiles(x) for x in list(cluster_df[1].CanonicalSMILES)]
            mcs1 = rdFMCS.FindMCS(mols, threshold=0.7)
            m1 = Chem.MolFromSmarts(mcs1.smartsString)
            img = Draw.MolToImage(m1, legend=cluster_df[0])
            img.show()



    def plot_experimental_parameter_distribution(self, df, reported_value_name):
        """
        This function first draws a pairplot of all environmenmental parameters specified in
        DataStructure.experimental_parameter_names. Second, it provides a distribution plot for each specified
        parameter separtatly

        @param df: DataFrame containing the environmental parameters and the reported endpoint
        @param reported_value_name: reported endpoint (e.g., reported logDT50 of individual experiments)
        """
        assert 'DataStructure' in self.object_name, \
            "This function cannot be applied to the object {}".format(self.object_name)
        sns.set(rc={"figure.figsize": (15, 15)})
        sns.set_theme(style="whitegrid")
        sns.set_style("ticks")
        sns.set_context(self.context)

        # draw pairplot for general overview and outlier detection
        print('Drawing pairplot...')
        output_file_path = os.path.join(self.get_data_directory(),
                                 'pairplot_parameter_distribution_{}_{}_{}.pdf'.format(self.data_type, self.tag,
                                                                              self.setup_name))
        sns.pairplot(df, plot_kws={'s': 20}) # hue can be added, but only works with categories
        print(f'Saving figure to {output_file_path}')
        plt.savefig(output_file_path)
        plt.close()

        # Draw distribution (violinplot + scatter) for each parameter, indicating it's unit(s) as hue
        print('Drawing violinplot of parameter distributions...')
        output_file_path = os.path.join(self.get_data_directory(),
                                 'violinplot_parameter_distribution_{}_{}_{}.pdf'.format(self.data_type, self.tag,
                                                                              self.setup_name))
        params = self.object.experimental_parameter_names
        sns.set(rc={"figure.figsize": (10, len(params))})
        sns.set_theme(style="whitegrid")
        sns.set_style("ticks")
        fig, axs = plt.subplots(int(round(len(params)/2)), 2)
        legend = True
        for ax, param in zip(fig.axes, params):
            print('\t- '+param)
            # todo: add units to names - needs to be checked first
            this_df = df.loc[:,[param, reported_value_name]]

            sns.stripplot(data=this_df, x=param, hue=reported_value_name, ax=ax, jitter=True,
                          palette="icefire", alpha=0.7, size=2, dodge=True, zorder=1, legend=legend)
            if legend:
                ax.legend(ncol=6, markerscale=3, title='logDT50 [log(days)]:', frameon=False, loc="lower left",
                          bbox_to_anchor=(0, 1.1), alignment='left', handletextpad=0.1, labelspacing=0.8,
                          borderpad=0, handlelength = 1.7)
                legend = False
            sns.violinplot(data=this_df, x=param, ax=ax, color='grey', legend=False, fill=False)

        print(f'Saving figure to {output_file_path}')
        plt.tight_layout()
        plt.savefig(output_file_path)
        plt.close()
