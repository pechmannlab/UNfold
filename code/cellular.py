import os, sys
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import networkx as nx
import subprocess
import glob
from Bio import AlignIO
import pickle


from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score


#--------------------------------------------------------------------
def QCDpred(seq):
    """
    Rund QCDpred to predict S.cerevisiae degrons
    This algorithm is developped by and taken from Johansson et al.:
    [Johansson et al. JMB (2023) doi.org/10.1016/j.jmb.2022.167915]
    """
    # Logistic regression model with L2 regularization (lambda=0.001)
    # Trained on WT data with >50 raw counts total in 4 the bins
    # Peptides classified degrons if PSI<2.2 and stable if PSI>2.8
    # Parameter are per amino acid counts in peptides of length 17
    model = {}
    model["intersect"] = -0.89102423
    model["C"] =   0.37721431
    model["D"] =  -0.78986558
    model["E"] =  -0.65124014
    model["K"] =  -0.15518666
    model["R"] =  -0.02030300
    model["H"] =  -0.02110156
    model["N"] =  -0.32782161
    model["Q"] =  -0.17676485
    model["A"] =   0.10844211
    model["G"] =  -0.37594135
    model["S"] =  -0.09627044
    model["T"] =  -0.08533912
    model["V"] =   0.43746326
    model["M"] =   0.31182498
    model["L"] =   0.53427787
    model["I"] =   0.61465146
    model["F"] =   0.52882600
    model["Y"] =   0.45253658
    model["W"] =   0.58693535
    model["P"] =  -0.25880796


    def calc_degron_prob(seq_list):
        p_list = []
        for seq in seq_list:
            # This model is only valid for peptides of length 17
            assert len(seq) == 17
            lin = model["intersect"] + sum([model[aa] for aa in seq])
            responce = 1/(1+np.exp(-lin))
            p_list.append(responce)
        return(p_list)



    def tile_sequence(seq, tile_len, tile_offset=1, add_ct_tile=True):
        """Cut sequence into tiles of length tile_len"""
        nres = len(seq)
        tile_begin = 0
        tile_end = tile_len
        tile_list = []
        while(tile_end <= nres):
            tile_list.append(seq[tile_begin:tile_end])
            tile_begin = tile_begin + tile_offset
            tile_end = tile_begin + tile_len
        if tile_end !=  nres+1 and add_ct_tile:
            tile_list.append(seq[-tile_len:])
        return(tile_list)

    tile_seq = tile_sequence(seq, 17, 1)
    tile_degron_prob = calc_degron_prob(tile_seq)
   
    resi = 9
    qcdDF = pd.DataFrame(columns=['tile_seq', 'tile_degron_prob', 'tile_aa', 'resi'])
    for i in range(len(tile_seq)):
        qcdDF.loc[(len(qcdDF))] = (tile_seq[i], tile_degron_prob[i], tile_seq[i][8], resi)
        resi += 1
   
    return qcdDF






if __name__ == '__main__':



    pdblist = [pdb.split('/')[-1].split('.')[0] for pdb in glob.glob('../data/pdb/*.pdb')]
    print(pdblist)



    # load data + compute degron stats
    mapDF = pd.read_csv("../data/processed/c3718_relpositions.txt", header=0, index_col=None, sep='\t')
    aln = pd.read_csv("../data/processed/c3718_aln.txt", header=None, index_col=0, sep='\t') 
    ss = pd.read_csv("../data/processed/consensus_aln_ss.txt", header=0, index_col=None, sep='\t')
    bpinit = pd.read_csv("../data/processed/assign_bpinit.txt", header=0, index_col=None, sep='\t')


    orflist = pd.read_csv("../data/processed/orflist.txt", header=0, index_col=None, sep='\t')
    substr = orflist[['orf', 'substrate']]
    tmp_substr = np.array(substr['substrate'])
    tmp_substr[tmp_substr=="+"] = 1
    tmp_substr[tmp_substr=="-"] = 0
    substr['substrate'] = list(tmp_substr)
   

    N_aln = list(ss[ss['element']=='b1']['aln'])[0]		# Nter before aln position of b1
    C_aln = list(ss[ss['element']=='a5']['aln'])[-1]	# Cter after aln position of a5
    terex = 10  # a4 plus loopy until b6

    degronDF = pd.DataFrame(columns=['pdb', 'substr', 'N', 'C', 'Nunfold', 'Cunfold', 'mean', 'mean_N', 'mean_C', 'mean_Np', 'mean_Cp', 'd70_N', 'd70_C', 'd70_Np', 'd70_Cp'])
    for px, pdb in enumerate(pdblist):

        current_substr = list(substr[substr['orf']==pdb]['substrate'])[0]  

        current_N = mapDF[mapDF['aln']==N_aln][pdb].values[0]
        current_C = mapDF[mapDF['aln']==C_aln][pdb].values[0]

        current_seq = np.array(aln.loc[pdb])
        current_aa = "".join(list( current_seq[current_seq!="-"] ))
    
        current_qcd = np.zeros(( len(current_aa) )) * np.nan      
        current_qcdDF = QCDpred(current_aa)
        current_start = np.array(current_qcdDF['resi'])[0] - 1
        current_qcd[ current_start:(current_start+int(len(current_qcdDF)) )  ] = np.array(current_qcdDF['tile_degron_prob'])
        np.savetxt("../data/SI/degron_"+pdb+".txt", current_qcd, fmt='%.3f')

        current_mean = np.around(np.nanmean(current_qcd), 3)
        current_mean_N = np.around(np.nanmean( current_qcd[:(current_N+terex)] ), 3)   
        current_mean_C = np.around(np.nanmean( current_qcd[(current_C-terex-1):current_C] ), 3) 

        current_mean_Np = np.around(np.nanmean( current_qcd[current_N:(current_N+terex)] ), 3)   
        current_mean_Cp = np.around(np.nanmean( current_qcd[(current_C-terex-1):] ), 3) 

        #50% 0.28652 # quantiles of all yeast tiles, data from KLL paper
        #70% 0.48584 
        # average degron probability: average of predicted scores
        # degron score: sum of positions with degron prob > 0.48584 (> 70% degron prob)

        d70_N = current_qcd[current_N:(current_N+terex)]
        current_d70_N = np.nansum( d70_N[d70_N > 0.48584] )

        d70_C = current_qcd[(current_C-terex-1):current_C]
        current_d70_C = np.nansum( d70_C[d70_C > 0.48584] )

        d70_Np = current_qcd[:(current_N+terex)]
        current_d70_Np = np.nansum( d70_Np[d70_Np > 0.48584] )

        d70_Cp = current_qcd[(current_C-terex-1):]
        current_d70_Cp = np.nansum( d70_Cp[d70_Cp > 0.48584] )


        current_bpinit = bpinit[bpinit['protein']==pdb]
        current_a2 = 0
        current_a5 = 0
            
        current_a2 += np.sum( np.array(current_bpinit['class']) == "a1" )
        current_a2 += np.sum( np.array(current_bpinit['class']) == "a2/b2" )
        current_a5 += np.sum( np.array(current_bpinit['class']) == "a3" )
        current_a5 += np.sum( np.array(current_bpinit['class']) == "a4/b5" )
        current_a5 += np.sum( np.array(current_bpinit['class']) == "a5" )

        current_a2 /= len(list(current_bpinit['class']))
        current_a5 /= len(list(current_bpinit['class']))

        Nter = 37  	# alignment position before start of b1
        Cter = 391	# alignment position after end of a5
        current_aln = np.array( aln.loc[pdb] )
        N = np.sum(current_aln[:Nter] != '-')
        C = np.sum(current_aln[Cter:] != '-')
 
        degronDF.loc[len(degronDF)] = (pdb, current_substr, N, C, np.around(current_a2, 3), np.around(current_a5,3), current_mean, current_mean_N, current_mean_C, current_mean_Np, current_mean_Cp, current_d70_N, current_d70_C, current_d70_Np, current_d70_Cp)

    degronDF.to_csv("../data/processed/degron_probs.txt", header=True, index=False, sep='\t')
    print(degronDF)


    # load protein abundance data
    paxdb = pd.read_csv("../data/yeast/yeast_paxdb.txt", header=11, index_col=None, sep='\t')
    listID = []
    for i in list(paxdb['string_external_id']):
        new_id = i.replace('4932.', '')
        listID.append(new_id)
    paxdb['orf'] = listID

    # load data from previous results & combine in data frame as features
    surf = pd.read_csv("../data/processed/surfaceAP_trajectory.txt", header=0, index_col=None, sep='\t')
    qf = pd.read_csv("../data/processed/qf_trajectory.txt", header=0, index_col=None, sep='\t')
    imdf = pd.read_csv("../data/processed/unfolding_intermediates.txt", header=0, index_col=None, sep='\t')
    Tf = pd.read_csv("../data/processed/Tf_contacts2.txt", header=0, index_col=None, sep='\t')

    gmxDF = pd.DataFrame(columns=['pdb', 'substr', 'abundance', 'agg', 'unfold', 'Tf', 'N_aln'])
    for px, pdb in enumerate(pdblist):
      
        current_pax = list(paxdb[paxdb['orf']==pdb]['abundance'])[0]
        current_substr = list(substr[substr['orf']==pdb]['substrate'])[0]

        current_traj = list(bpinit[bpinit['protein']==pdb]['traj'])
        list_agg = []
        list_intmi = []
        for i in current_traj:
            current_surf = np.array(surf[i])
            current_qf = np.array(qf[i])
            current_unfolded = current_qf < 0.3
            current_agg = current_surf[current_unfolded]
            list_agg += list(current_agg)

            current_intmi = imdf[imdf['traj']==i]['end'] - imdf[imdf['traj']==i]['start']
            #current_intmi = imdf[imdf['traj']==i]['mindist']
            list_intmi.append(int(current_intmi))
        current_meanagg = np.mean( np.array(list_agg) )
        #current_intmi = np.mean( np.array(list_intmi))
        current_intmi = np.sum(np.array(list_intmi) > 200)/len(list_intmi)
        current_Tf = list(Tf[Tf['pdb']==pdb]['Tf'])[0]
        current_Naln = list(Tf[Tf['pdb']==pdb]['N_aln'])[0]

        gmxDF.loc[len(gmxDF)] = (pdb, current_substr, current_pax, current_meanagg, current_intmi, current_Tf, current_Naln)

    print(gmxDF)
    gmxDF.to_csv("../data/processed/genomix.txt", header=True, index=False, sep='\t')





    ### UNFOLD PREDS
    data = gmxDF[['substr', 'abundance', 'agg', 'unfold', 'Tf', 'N_aln']]
    data['abundance'] = np.log(data['abundance'])
    
    # 1. only fullfit/predict
    y = np.array(data['unfold'], dtype=float)
    X = np.array(data[['abundance', 'agg', 'substr', 'Tf', 'N_aln']])

    regr = RandomForestRegressor(max_depth=2, random_state=0)
    regr.fit(X, y)
    pred = regr.predict( X )
    resDF = pd.DataFrame({'true': list(y), 'pred': list(pred) })
    resDF.to_csv("../data/processed/unfold_fullfit.txt", header=True, index=False, sep='\t')
    np.savetxt("../data/processed/unfold_fullfit_fa.txt", regr.feature_importances_)

    # 2. leave-one-out cross-validation for predictive power
    # !! omit one protein from training and discuss !! too few data points for proper training ... 
    X0 = np.array(data.iloc[0][['abundance', 'agg', 'substr', 'Tf', 'N_aln']]).reshape(1, -1)
    y0 = np.array(data['unfold'])[0]
    data = data.drop([0], axis=0)
 
    y = np.array(data['unfold'], dtype=float)
    X = np.array(data[['abundance', 'agg', 'substr', 'Tf', 'N_aln']])

    pred = np.zeros(( len(y) ))
    pred0 = np.zeros(( len(y) ))
    fa = np.zeros(( len(y), np.shape(X)[1] ))
    for i in range( len(y) ):          # len(pdblist)
        sel_train = np.ones(( len(y) ), dtype=bool)
        sel_train[i] = False

        y_train = y[sel_train]
        X_train = X[sel_train,:]
        X_pred = X[i,:].reshape(1, -1)

        regr = RandomForestRegressor(max_depth=2, random_state=0)
        regr.fit(X_train, y_train)
        pred[i] = regr.predict( X_pred )
        pred0[i] = regr.predict( X0 )

        fa[i,:] = regr.feature_importances_

    p0 = np.mean(pred0)

    y_true = np.zeros(( 16 ))
    y_true[0] = y0
    y_true[1:] = y

    y_pred = np.zeros(( 16 ))
    y_pred[0] = p0
    y_pred[1:] = pred

    resDF = pd.DataFrame({'true': list(y_true), 'pred': list(y_pred) })
    resDF.to_csv("../data/processed/unfold_xval.txt", header=True, index=False, sep='\t')
    np.savetxt("../data/processed/unfold_xval_fa.txt", fa)





    ### SUBSTR PREDS
    # 1. only fullfit/predict
    data = gmxDF[['substr', 'abundance', 'agg', 'unfold', 'N_aln']]
    data['abundance'] = np.log(data['abundance'])
 
    y = np.array(data['substr'], dtype=int)
    X = np.array(data[['abundance', 'agg', 'unfold', 'N_aln']])

    rfc = RandomForestClassifier(n_estimators=100, random_state = 0, n_jobs = -1)
    rfc.fit(X, y)
    pred = rfc.predict( X )
    resDF = pd.DataFrame({'true': list(y), 'pred': list(pred) })
    resDF.to_csv("../data/processed/substr_fullfit.txt", header=True, index=False, sep='\t')
    np.savetxt("../data/processed/substr_fullfit_fa.txt", rfc.feature_importances_)


    # 2. leave-one-out cross-validation for predictive power
    pred = np.zeros(( len(y) ))
    pred0 = np.zeros(( len(y) ))
    fa = np.zeros(( len(y), np.shape(X)[1] ))
    for i in range( len(y) ):          # len(pdblist)
        sel_train = np.ones(( len(y) ), dtype=bool)
        sel_train[i] = False

        y_train = y[sel_train]
        X_train = X[sel_train,:]
        X_pred = X[i,:].reshape(1, -1)

        rfc = RandomForestClassifier(n_estimators=100, random_state = 0, n_jobs = -1)
        rfc.fit(X_train, y_train)
        pred[i] = rfc.predict( X_pred )

        fa[i,:] = rfc.feature_importances_

    resDF = pd.DataFrame({'true': list(y), 'pred': list(pred) })
    resDF.to_csv("../data/processed/substr_xval.txt", header=True, index=False, sep='\t')
    np.savetxt("../data/processed/substr_xval_fa.txt", fa)

    ac = accuracy_score(y, pred)
    print(ac)
