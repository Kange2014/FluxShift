import sys, os

from IPython import get_ipython
from IPython.display import clear_output
get_ipython().run_line_magic('matplotlib','inline')

from matplotlib import pyplot as plt

# Panda python module for dataframe and data storage/manipulation
import pandas as pd
pd.set_option('mode.use_inf_as_null',True)
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 999)
pd.set_option('precision', 3)

import seaborn as sns
sns.set(style="white")
c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 = sns.color_palette("Set2", 10)
b1, b2, b3, b4, b5, b6 = sns.color_palette("Blues")

from copy import deepcopy
import numpy as np
import scipy.io
import scipy.sparse as sparse
import random as rand

# For MCMC sampling analysis

def remove_unwanted_rxns(df,reactions_to_remove):
    df_reduced = df.copy(deep=True)
    
    for i in df.columns:
        if i in reactions_to_remove:
            df_reduced.drop(str(i), axis=1, inplace=True)
            
    return df_reduced

def prune_mcmc_fluxes(df,m):
    
    #1) normalize to glc uptake rate
    
    for react in df.columns:
        for i in df[react].index[0:1]:
            if i == 0:
                glc_uptake = df.EX_glc__D_e[i]
            
            df[react][i] = np.true_divide(df[react][i], glc_uptake)
            
    #2) weed out reactions with flux over 10*glc uptake rate for over 75% of the flux states
    # produces two new df with removed columns corresponding to list rxn_ex_rpp_removed and
    # those with fluxes > 10 from 1-25% of the time
    
    rxn_removed = []
    read_to_struct = []
    read_to_struct2 = []
    
    for react in df.columns:
        if len(df[abs(df[react])<10]) < len(df.index)*.75:
            rxn_removed.append(react)
            read_to_struct.append({ 'react_id':react, 'subsystem':m.reactions.get_by_id(react).subsystem,  'percent of states with flux < 10':len(df[abs(df[react])<10]), 'min abs flux':abs(df[react]).min(), 'max abs flux':abs(df[react]).max()})

        # see how many remaining reactions have flux > 10 for 1-25% of the flux distributions
        elif len(df[abs(df[react])<10]) != len(df.index):
            read_to_struct2.append({'react_id':react, 'subsystem':m.reactions.get_by_id(react).subsystem,  'percent of states with flux < 10':len(df[abs(df[react])<10]), 'min abs flux':abs(df[react]).min(), 'max abs flux':abs(df[react]).max()})


    df_rxns_high_flux_removed = pd.DataFrame(read_to_struct)  
    df_rxns_high_flux_partial = pd.DataFrame(read_to_struct2)  

    df = df[list(set(df.columns) - set(rxn_removed))]
    
    return (df, df_rxns_high_flux_removed, df_rxns_high_flux_partial)


def get_mean_std_pvalue(sample1, sample2, m, state1, state2):    
    ''' takes in two dataframes that are sampled flux distributions (sample1 and sample2 are two sampled states)
    m is an M model and state1 and state2 are string names of the sampled states. It outputs a new dataframe with
    statistical measures such as means, medians and std's for each state and the pvalue which can be used to determine
    how much overlap there is between two flux distributions of a reaction in each respective state'''
    
    out = []
    
    for r in m.reactions:
        r = r.id
        if r in sample1.columns:
            cond1 = sample1[r]
        else:
            cond1 = pd.Series(np.zeros(10000))
        if r in sample2.columns:
            cond2 = sample2[r]
        else:
            cond2 = pd.Series(np.zeros(10000))


        vals = []
        for i in range(0,10):
            cond1 = cond1.reindex(np.random.permutation(cond1.index))
            cond2 = cond2.reindex(np.random.permutation(cond2.index))
            z = cond1 - cond2
            x = len(z[z>0]) + 1
            y = len(z[z<0]) + 1
            k = min(x,y)
            vals.append(k)
        p = np.mean(vals)/len(cond1)*2 ## is it vals or k?
        c1mean = cond1.mean()
        c2mean = cond2.mean()
        c1median = cond1.median()
        c2median = cond2.median()
        c1std = cond1.std()
        c2std = cond2.std()
        
        out.append({'reaction':r,'p%dmean'%state1:c1mean, 'p%dmean'%state2:c2mean,'p%dmedian'%state1:c1median, 'p%dmedian'%state2:c2median,'p%dstd'%state1:c1std,'p%dstd'%state2:c2std,'pval':p})
        
    out = pd.DataFrame(out)
    return out

def plot_rxn_shifts_all_strains(df, df2, r, xmax, xmin, ymax):
    
    name_graph = r+'_'+'all_strains'+'_rxn_shifts'
    
    fig, axes = plt.subplots(nrows=1, ncols=2)  # change for adding in bis
    
    #c1, c2, c3 = sns.color_palette("Set1", 3)
    
    c1, c2 = sns.color_palette("Set2", 2)

    keys = df.keys()
    print keys
    labels = []
    if keys[0] == 'Control':
        labels.append('phase 1')
        labels.append('phase 2')
    else:
        labels.append('Control')
        labels.append('20aas')
        
    df_tmp = df[keys[0]]
    df_tmp2 = df2[keys[0]]
    strain = keys[0]
    
    print df_tmp.shape
    
    #df_tmp[r].plot(kind='density',label=1, ax=axes[0,0],color=sns.color_palette()[0], title=strain)
    #df_tmp2[r].plot(kind='density',label=2,ax=axes[0,0],color=sns.color_palette()[1], title=strain)
    sns.kdeplot(df_tmp[r], shade=True, color=c2, alpha=0.7,ax=axes[0],label= labels[0]);
    sns.kdeplot(df_tmp2[r], shade=True, color=c3, alpha=0.7,ax=axes[0],label= labels[1])
    axes[0].set_xlim(xmin,xmax)
    axes[0].set_ylim(0,ymax)
    axes[0].set_title(strain)
    
    df_tmp = df[keys[1]]
    df_tmp2 = df2[keys[1]]
    strain = keys[1]
    
    sns.kdeplot(df_tmp[r], shade=True, color=c2, alpha=0.7,ax=axes[1],label= labels[0]);
    sns.kdeplot(df_tmp2[r], shade=True, color=c3, alpha=0.7,ax=axes[1],label= labels[1])
    axes[1].set_xlim(xmin,xmax)
    axes[1].set_ylim(0,ymax)
    axes[1].set_title(strain)
    
    fig.tight_layout()

    fig = plt.gcf()
    fig.set_figwidth(17)
    fig.set_figheight(10)
	
# Classify flux shifts

import math

def gen_rxn_to_subsystem(m):
    out = []
    for r in m.reactions:
            if 'sink' in r.name.lower() or 'exchange' in r.name.lower():
                subsystem = "exchange/sink reaction"
            elif len(r.subsystem) == 0:
                subsystem = "None"
            else:
                subsystem = r.subsystem
            out.append({'reaction':r.id, 'subsystem':subsystem})
    return pd.DataFrame(out)

def gen_rxn_to_genes(m):
    out = []
    for r in m.reactions:
        if len(r._genes) <1:
                out.append({'reaction':r.id, 'gene':np.nan})
        else:
            for g in r._genes:
                out.append({'reaction':r.id, 'gene':g.id})
            
    return pd.DataFrame(out)

def classify_flux_change(x):
    flux_change_pval = 0.05 
    abs_mean_rel_diff = 0.105 # 5% of data  
    abs_kde_rel_diff = 0.0850 # 5% of data
    
    if x['pval'] < flux_change_pval and ( x['abs_mean_rel_diff'] > abs_mean_rel_diff):
        return 'change'
    else:   
        return 'same'

def find_flux_shifts(df,m,state1,state2):
    print 'no. rxns before filtering out exchange rxns:', len(df.index)

    r2g = gen_rxn_to_genes(m)
    r2s = gen_rxn_to_subsystem(m)

    # remove certain subsystems
    subsystems_to_omit = ['Transport, Inner Membrane','Transport, Outer Membrane Porin','Transport, Outer Membrane']
    df = pd.merge(df, r2s, on='reaction')
    df = df[~df.subsystem.isin(subsystems_to_omit)]

    print "no. rxns after filtering out exchange rxns:", len(df.index)

    # round fluxes to 5 decimal places
    dec = 5    
    if state1 == 1:
        df['p1mean'] = np.round(df['p1mean'], decimals=dec) 
        df['p2mean'] = np.round(df['p2mean'], decimals=dec)
      
    else:
        df['p2mean'] = np.round(df['p2mean'], decimals=dec) 
        df['p3mean'] = np.round(df['p3mean'], decimals=dec)

    
    # absolute fold difference between means of states
    if state1 == 1:
        abs_mean_fold_diff = df['p2mean']/df['p1mean']         # fold difference (to do: should be absolute???)
    else:
        abs_mean_fold_diff = df['p3mean']/df['p2mean']         # fold difference (to do: should be absolute???)
        
    abs_mean_fold_diff[abs_mean_fold_diff==np.inf] = 1000      # change infinite vs zero's so they can be dealt with
    abs_mean_fold_diff[abs_mean_fold_diff==0] = 0.001
    abs_mean_fold_diff = abs_mean_fold_diff.fillna(1) # fill NaN with 1's for fold changes
    df['abs_mean_fold_diff'] = abs_mean_fold_diff
    #abs_mean_fold_diff.sort(ascending=False) 
    
    # absolute relative mean differences
    if state1 == 1:
        abs_mean_rel_diff = abs(abs(df['p2mean'])-abs(df['p1mean'])) # absolute difference in flux change
    else:
        abs_mean_rel_diff = abs(abs(df['p3mean'])-abs(df['p2mean'])) # absolute difference in flux change
    df['abs_mean_rel_diff'] = abs_mean_rel_diff
            
    # classify shifts
    df['flux_change'] = df.apply(lambda x: classify_flux_change(x), axis=1)
    
    # create separate dataframes for fluxes that stay the same vs shift
    change = df[df['flux_change'] == 'change']
    same = df[df['flux_change'] == 'same']

    col_to_omit =[]
    for i in change['reaction'].index:
        if 'EX' in change['reaction'][i] or 'pp' in change['reaction'][i]:
            col_to_omit.append(change['reaction'][i])
        
    same = same[~same.reaction.isin(col_to_omit)]
    change = change[~change.reaction.isin(col_to_omit)]

    print '\nflux change counts:'
    print 'dont shift', len(same)
    print 'do shift', len(change), '\n'
    
    return (df,same,change)

	
def get_zscore(sample1, sample2, m, flux_shift_stats_dict, strain):    
    ''' takes in two dataframes that are sampled flux distributions (sample1 and sample2 are two sampled states)
    m is an M model. It outputs a new dataframe with reaction z scores'''
    
    out = []
    counter = 0
    
    # compute z score per reaction,x, for two conditions (p = population)
    # z = x- mean(p) / std(p)
    
    # get population mean difference, mu_p (mean(p))
    mu_p = flux_shift_stats_dict[strain]['abs_mean_rel_diff'].mean()
    std_p = flux_shift_stats_dict[strain]['abs_mean_rel_diff'].std()
    
    # get x:
    for r_j in m.reactions:
        
        clear_output(wait=True)
        print('***PROGRESS: %d/%d reactions***\n' % (counter, len(m.reactions)))
        counter += 1
        sys.stdout.flush()        
        
        z_score_Rj = [] # z_score_Rj is a list of all z_scores for all reactions, r_j in N_j
                
        for i in range(0,10):
            react = r_j.id
                
            if react in sample1.columns:
                cond1 = sample1[react]
            else:
                cond1 = pd.Series(np.zeros(10000))
            if react in sample2.columns:
                cond2 = sample2[react]
            else:
                cond2 = pd.Series(np.zeros(10000))

            # 10,000 flux differences randomly permuted!    
            cond1 = cond1.reindex(np.random.permutation(cond1.index))
            cond2 = cond2.reindex(np.random.permutation(cond2.index))
            
            z = abs(cond1 - cond2)
            denom = np.true_divide(z.std(),np.sqrt(len(z)))
            z_score_i = (np.true_divide(z.mean()-mu_p,std_p)) # this is reaction z score per reaction
            z_score_Rj.append(abs(z_score_i))
        
        z_score=np.mean(z_score_Rj) 
        
        out.append({'reaction':react,'z_score':z_score})
        
    out = pd.DataFrame(out)
    return out

def get_met_zscore(sample1, sample2, m, flux_shift_stats_dict, DF_reaction_zscore, strain):
    read_to_struct = []
    m_met_Nj =[]
    counter = 0
    
    # compute m_met,j for all metabolites from DF reaction zscores to get mu_p and std_p

    for j in m.metabolites:
        m_met = 0
        N_j = len(m.metabolites.get_by_id(j.id).get_reaction()) # N_j = number of reactions met_j is involved in
        if N_j != 0:  
            for r_j in m.metabolites.get_by_id(j.id).get_reaction():
                react = r_j.id
                m_met += DF_reaction_zscore[strain].z_score[DF_reaction_zscore[strain].reaction == react].values[0]
                
            m_met_Nj.append(np.true_divide(m_met,np.sqrt(N_j)))
            
    mu_p = np.mean(m_met_Nj)
    std_p = np.std(m_met_Nj)
    #print mu_p, std_p, len(m_met_Nj)

    mu_r_p = flux_shift_stats_dict[strain]['abs_mean_rel_diff'].mean()
    std_r_p = flux_shift_stats_dict[strain]['abs_mean_rel_diff'].std()
                
    # compute m_met,j for 1000 randomly generated sampled reaction zscores

    for j in m.metabolites:    
        m_met_Nj=[]  # m_met_Nj is a list that holds all m_met_j z_scores from 100,000 randomly generated flux difference samples
        N_j = len(m.metabolites.get_by_id(j.id).get_reaction()) # N_j = number of reactions met_j is involved in
    
        if N_j != 0:    
            clear_output(wait=True)
            print('***PROGRESS: %d/%d metabolites***\n' % (counter, len(m.metabolites)))
            counter += 1
            sys.stdout.flush()
            
            for i in range(0,10):
                
                z_score_Rj = [] # z_score_Rj is a list of all z_scores for all reactions, r_j in N_j
                
                for r_j in m.metabolites.get_by_id(j.id).get_reaction():
                    react = r_j.id
                
                    if react in sample1.columns:
                        cond1 = sample1[react]
                    else:
                        cond1 = pd.Series(np.zeros(10000))
                    if react in sample2.columns:
                        cond2 = sample2[react]
                    else:
                        cond2 = pd.Series(np.zeros(10000))

                    # 10,000 flux differences randomly permuted!    
                    cond1 = cond1.reindex(np.random.permutation(cond1.index))
                    cond2 = cond2.reindex(np.random.permutation(cond2.index))
            
                    z = abs(cond1 - cond2)
                    denom = np.true_divide(z.std(),np.sqrt(len(z)))
                    z_score_i = (np.true_divide(abs(z.mean()-mu_r_p),std_r_p)) # this is my reaction z score per reaction
                    z_score_Rj.append(z_score_i)
                    
                m_met_Nj.append(np.true_divide(np.sum(z_score_Rj),np.sqrt(N_j)))
                        
        mean_m_met_Nj = np.mean(m_met_Nj)        
                
        Z_met_j = np.true_divide(abs(mean_m_met_Nj - mu_p),std_p)
        #print j, m_met, mean_m_met_Nj, std_m_met_Nj, mu_met_Nj #, np.mean(m_met_Nj), np.std(m_met_Nj)
        read_to_struct.append({'metabolite':j.id,'num_rxns_involved':N_j, 'zscore':Z_met_j})
            
    return pd.DataFrame(read_to_struct)
