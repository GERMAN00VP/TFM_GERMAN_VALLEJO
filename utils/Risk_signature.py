import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

from sksurv.linear_model import CoxnetSurvivalAnalysis
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter

def plot_coefficients(coefs, n_highlight=5,save=False,show=False):
    """

    Args:
        coefs (pd.dataframe): dataframe with the lasso coefficients of each gene along alphs values
        n_highlight (int): number of genes highlighted
        save (bool|str, optional): path where the plot is saved if save is a str. Defaults to False.
        show (bool, optional): if True it shows the plot. Defaults to False.
                
    """
    _, ax = plt.subplots(figsize=(9, 6))

    alphas = coefs.columns
    for row in coefs.itertuples():
        ax.semilogx(alphas, row[1:], ".-", label=row.Index)

    alpha_min = alphas.min()
    top_coefs = coefs.loc[:, alpha_min].map(abs).sort_values().tail(n_highlight)
    for name in top_coefs.index:
        coef = coefs.loc[name, alpha_min]
        plt.text(
            alpha_min, coef, name + "   ",
            horizontalalignment="right",
            verticalalignment="center"
        )

    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.grid(False)
    ax.set_xlabel("Alpha")
    ax.set_ylabel("Coefficient")
    
    if save:
        plt.savefig(f"{save}alpha_coefficients_plot.png",bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()
    

def plot_alpha_selection(gcv,save=False,show=False):
    
    """
    
    Args:
        gcv (sklearn.model_selection._search.GridSearchCV): the result of a grid search, used for extracting all the values nedeed for the plot.
        save (bool|str, optional): path where the plot is saved if save is a str. Defaults to False.
        show (bool, optional): if True it shows the plot. Defaults to False.
        
    """
    
    cv_results = pd.DataFrame(gcv.cv_results_) 
    alphas = cv_results.param_coxnetsurvivalanalysis__alphas.map(lambda x: x[0])
    mean = cv_results.mean_test_score
    std = cv_results.std_test_score

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(alphas, mean)
    ax.fill_between(alphas, mean - std, mean + std, alpha=.15)
    ax.set_xscale("log")
    ax.set_ylabel("concordance index")
    ax.set_xlabel("alpha")
    ax.axvline(gcv.best_params_["coxnetsurvivalanalysis__alphas"][0], c="C1")
    ax.axhline(0.5, color="grey", linestyle="--")
    ax.grid(True)
    
    if save:
        plt.savefig(f"{save}alpha_selection_plot.png",bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()
        
    
def plot_non0_coefs(gcv,genes,save=False,show=False):
    """

    Args:
        gcv (sklearn.model_selection._search.GridSearchCV): the result of a grid search, used for extracting all the values nedeed for the plot.
        genes (list): list of genes (with the same format that those in training set columns)
        save (bool|str, optional): path where the plot is saved if save is a str. Defaults to False.
        show (bool, optional): if True it shows the plot. Defaults to False.
    
    Returns:
        non_zero_coefs (pandas.DataFrame): df with genes selected by the lasso cox in the best estimator of the grid search
    """
    
    best_model = gcv.best_estimator_.named_steps["coxnetsurvivalanalysis"]
    best_coefs = pd.DataFrame(
        best_model.coef_,
        index=genes,
        columns=["coefficient"]
    )

    non_zero = np.sum(best_coefs.iloc[:, 0] != 0)
    print(f"Number of non-zero coefficients: {non_zero}")

    non_zero_coefs = best_coefs.query("coefficient != 0")
    coef_order = non_zero_coefs.abs().sort_values("coefficient").index

    _, ax = plt.subplots(figsize=(6, 8))
    non_zero_coefs.loc[coef_order].plot.barh(ax=ax, legend=False)
    ax.set_xlabel("coefficient")
    ax.grid(True)
    
    if save:
        plt.savefig(f"{save}non0_coefs_plot.png",bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()
        
    return non_zero_coefs

def Risk_Sign_Calculator(df, LASSOgenes, LASSOcoefs,set_name="Training",save=False,show=False):

    """
    
    Args: 

        df (pd.DataFrame): (samples, genes + [OS,OS.time]).
        LASSOgenes (list|np.array|pd.Series): genes in the risk signature.
        LASSOgenes (list|np.array|pd.Series): coefficients for the genes in the risk signature.
        set_name (str): name of the dataset passed. Defaults to 'Training'. 
        save (bool|str, optional): path where the plot is saved if save is a str. Defaults to False.
        show (bool, optional): if True it shows the plot. Defaults to False.
    
    Description: 
        Calculates the value of the risk for each sample, stratificates the patients in low and hight risk and plot the Kaplan-Meier curve.
    """
    
    print(f"Calculating {set_name} set risk score")
    
    for i in range(len(LASSOgenes)):
        df[f'{LASSOgenes[i]}_coef'] = df[f'{LASSOgenes[i]}']*LASSOcoefs[i]
    
    df['Risk_score'] = df[[gene_name + '_coef' for gene_name in LASSOgenes]].sum(axis = 1)
    
    median = df['Risk_score'].median()

    df['Risk_group']= ''
    df.loc[df['Risk_score'] > median, 'Risk_group'] = 1 # Esto no es siempre solo por la distribucion de los datos (normalmente es >median)
    df.loc[df['Risk_score'] < median, 'Risk_group'] = 0


    print(df.head())

    print()
    print(f'Risk Score Median: {median}')

    print(df.isnull().sum())

    df.dropna(inplace = True)
    print('Null values dropped.')

    ################################################################################################

    fig, ax = plt.subplots(1, 1, figsize = (6, 3))
    kmf = KaplanMeierFitter()
    T = df['OS.time'].loc[df['Risk_group'] == 1]
    E = df['OS'].loc[df['Risk_group'] == 1]
    kmf.fit(durations = T, event_observed = E, label = 'High')
    kmf.plot_survival_function()

    kmf1 = KaplanMeierFitter()
    T1 = df['OS.time'].loc[df['Risk_group'] == 0]
    E1 = df['OS'].loc[df['Risk_group'] == 0]
    kmf1.fit(durations = T1, event_observed = E1, label = 'Low')
    kmf1.plot_survival_function()

    results = logrank_test(T1,T,event_observed_A=E1, event_observed_B=E)

    add_at_risk_counts(kmf,kmf1, ax=ax, fig=fig, labels=['Low', 'High'], rows_to_show=['At risk'])

    fontsize = 13

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        #tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        #tick.label1.set_fontweight('bold')

    ax.text(0.2, 0.1, 'p.value= ' + str("{0:.6f}".format(results.p_value)), fontsize=10,weight='bold')
    fig.patch.set_facecolor('w')
    ax.spines[['right', 'top']].set_visible(False)
    plt.xlabel('Days passed')
    plt.title('Survival of different RS group')
    
    if save:
        plt.savefig(f"{save}KM_plot_{set_name}_set.png",bbox_inches="tight")
        
    if show:
        plt.show()


def create_results_folder(res_dirname):
    if not os.path.exists(res_dirname):
        os.mkdir(res_dirname)
    else:
        print("Path exists")


def risk_signature_maker(training_set,imput_gene_list, validation_sets={},out_dir="./Results_Risk_signature/",variation_filter=None,show=True,n_jobs=4):
    """

    Args:
        training_set (pandas.DataFrame): the df where the risk signature is trained shape = (samples,genes+[OS,OS.time]).Deduplicated and without missing values. 
        imput_gene_list (list): list of genes (with the same format that those in training set columns)
        validation_sets (dict): Defaults to {}. A dict like: {"validation_dataset_name":df}, dataframes must have same format as training_set
        out_dir (str, optional): _description_. Defaults to "./Results_Risk_signature/". Where to place the results, if it doesnt exists it's created. If False, does'nt save the results.
        variation_filter (float, optional): _description_. Defaults to None. If a value is provided, filter out the genes under that variance threshold in the thraining set.
        show (bool, optional): _description_. Defaults to True. Whether to sow the plots or not.
        
        
    Returns:
        dchp1 (pandas.Dataframe): a dataframe with the significative (p<0.05) genes and coeficcients of the risk signature.
        genes_tested: (list) the gene list that is finally tested after the intersection with training and validation datasets.
        
        
    Description:
        Creates a risk signature from a training set and a gene list, it stores in the out_dir the results_df and plots of the coeficient shrinkage of the lasso, of the alpha value selected, 
        the lasso genes with non 0 coefficients, and the KM plots on the training and validation datasets.
         
    """
    
    if out_dir:
        create_results_folder(out_dir) # Create the results folder
        if out_dir[-1]!="/":
            out_dir=out_dir+"/" #  To accurately save the results
    
    # Separate the OS and OS.time columns
    
    os_cols = ['OS','OS.time']
    training_set_os = training_set[os_cols]
    training_set.drop(os_cols,axis=1,inplace=True)
    
    # Filter out the genes under the variance threshold
    if variation_filter != None:
        if training_set[training_set.columns[training_set.var()>variation_filter]].shape[1] !=0:
            print(f"Filtering out genes with var < {variation_filter}")
            training_set = training_set[training_set.columns[training_set.var()>variation_filter]]
        else: 
            print("Var filtering not applied because all genes were filtered out, set a lower threshold or check if your data is transformed (like log1p() data).")
    
    validation = False
    # Then we need to do the intersection of the genes with the three datasets and the list of genes 
    if len(validation_sets)>0:
        validation_names = [name for name in validation_sets.keys()]
        validation=True
        
    genes =list (set(imput_gene_list).intersection(training_set.columns)) # Do the intersection with the training set
    
    if validation:
        for dataset in validation_names:
            genes =list (set(genes).intersection(validation_sets[dataset].columns)) # Do the intersection with the validation dfs
    
    # Create a bool os df of the training set
    training_set_os_bool= training_set_os.copy()
    training_set_os_bool["OS"] = training_set_os_bool.OS==1 # 1 == Decessed
    
    # Convert the os df to array with the correct format for CoxnetSurvivalAnalysis
    os_time_arr = np.array([(os,time) for os, time in zip(training_set_os_bool.OS,training_set_os_bool["OS.time"])], dtype=[('Status', '?'), ('Survival_in_days', '<f8')])
    
    
    # Do a previous cox lasso analysis
    X = StandardScaler().fit_transform(training_set[genes])    
    coxnet = CoxnetSurvivalAnalysis(l1_ratio=1.0, alpha_min_ratio=0.01).fit(X,os_time_arr)
    coefficients_lasso = pd.DataFrame(
        coxnet.coef_,
        index=genes,
        columns=np.round(coxnet.alphas_, 5))

    
    
    estimated_alphas = coxnet.alphas_ # Use the alphas of the previous analisis as imput for the param grid
    cv = KFold(n_splits=5, shuffle=True, random_state=0) # Define the cross validation object
    
    # Perform the grid search on the training set
    gcv = GridSearchCV(
        make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=1)),
        param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in estimated_alphas]},
        cv=cv,
        error_score=0.5,
        n_jobs=n_jobs).fit(training_set[genes], os_time_arr) 
    
    plot_alpha_selection(gcv, save=out_dir,show=show) # Plot the alpha selection curve
    
    non_zero_coefs = plot_non0_coefs(gcv,genes, save=out_dir,show=show) # Plot the non 0 coefficients
    
    non_0_coefs_lasso = pd.concat([training_set_os,training_set[non_zero_coefs.index]],axis=1) # Concatenate the OS and filtered training_set dataframes
    
    cph = CoxPHFitter().fit(non_0_coefs_lasso,'OS.time', 'OS') # Perform the multivariant cox analysis
    dcph = cph.summary[['coef', 'exp(coef)', 'p']] # Extract the meaningfull info from the cox result summary
    dchp1 = dcph.loc[dcph['p']<0.05] # Filter out the nos significative genes
    
    if out_dir:
        dchp1.to_csv(f"{out_dir}Multivariant_Cox_significative_results.csv")
        
    df_risk_signature = non_0_coefs_lasso[dchp1.index.tolist()+["OS","OS.time"]] # Define the dataframe for the risk signature calculation
    LASSOgenes = dchp1.index.tolist()
    LASSOcoefs = dchp1.coef.tolist() 
    print(f"The risk signature has a legth of: {len(LASSOgenes)} genes")
    
    Risk_Sign_Calculator(df_risk_signature, LASSOgenes, LASSOcoefs,set_name="Training",
                         save=out_dir,show=False) # Plot the training dataset Kaplan-Meier curve
    
    if validation: # Perform the validation
        for name in validation_names:
            df_risk_signature = validation_sets[name][dchp1.index.tolist()+["OS","OS.time"]]
            Risk_Sign_Calculator(df_risk_signature, LASSOgenes, LASSOcoefs,set_name=name,
                            save=out_dir,show=False) # Plot the validation datasets Kaplan-Meier curve
            
            
    plot_coefficients(coefficients_lasso, save=out_dir,show=show,n_highlight=5) # Plot the alphas coefficients shrinkage
    
    
    return dchp1,genes