import os, sys
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import subprocess
import glob
from Bio import AlignIO
import pickle

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import sbmOpenMM
import mdtraj as md
from Bio.PDB.PDBParser import PDBParser

from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def project_asa(alnfile, D, fileOUT):
    """
    compute and project profiles from input dictionaty onto alignment
    alnfile: input alignment file
    D: input dictionary, values per alignemnt position
    fileOUT: filename for output 
    """

    alignment = AlignIO.read(open(alnfile), "clustal")

    nrow = len(alignment)
    ncol = alignment.get_alignment_length()
    alnmat = np.zeros(( nrow, ncol ))
    list_substr = []
    list_id = []
    list_ix = []

    for ix, record in enumerate(alignment):
        seq_aln = np.array(record.seq)
        seq_ref = "".join(list(seq_aln[seq_aln!='-']))
        current_record = record.id.split('_')[0]

        if "dssp" not in record.id and current_record not in list_id and current_record != "space":
            list_id.append(current_record)
            list_ix.append(ix)

            # get profile from dictionary D
            ident = record.id.split('_')[0]
            asa = D[ident]
            asa_aln = np.zeros(( len(seq_aln) )) * np.nan

            pos_aln = 0
            pos_ref = 0
            while (pos_aln < ncol):
                if seq_aln[pos_aln] == '-':
                    pos_aln += 1
                else:
                    asa_aln[pos_aln] = asa[pos_ref]
                    pos_aln += 1
                    pos_ref += 1

            alnmat[ix,:] = asa_aln
 
    alnmat = alnmat[np.array(list_ix), :]
    alnmat = np.around(alnmat, 4)

    alndf = pd.DataFrame(alnmat, index=list_id)
    alndf.to_csv(fileOUT, sep='\t', header=True, index=False, na_rep="NA")

    return alndf, alnmat


 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_frustration_data(pdblist):
    """
    wrapper function to load frustration data precomputed with Frustratormeter2.0 (REF)
    input paths adapted from program output
    """

    frustdct = {}

    for PDB in pdblist:
        current_frust = "../data/frustration/"+PDB+".done/FrustrationData/"+PDB+".pdb_configurational_5adens"

        current_data = pd.read_csv(current_frust, header=0, index_col=False, sep=' ')
        frustdct[PDB] = np.array(current_data['relHighlyFrustrated'])

    return frustdct



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def color_bfac(pdbIN, dict_bfac, pdbOUT):
    """
    Function to write provided values to bfac column of PDB for coloring
    Here for single chain PDBs, no accounting for chainID
    dict_fact: dictionary of new bfac values for each residue id
    pdbIN: input PDB
    pdbOUT: output PDB
    """

    with open(pdbIN) as pdb_file:
        pdb = pdb_file.readlines()

    default_value = 0
    new_pdb = []
    for line in pdb:
        if line[0:6] == "ATOM  ":
            chain_id = line[20:22].strip()		# not using, but modify here for multi-chain PDBs
            res_id = int( line[23:26].strip() )
            bfac = dict_bfac.get(res_id, default_value)

            new_line = line[:60] + ("{bfac:6.2f}".format(bfac=bfac) ) + line[66:]
            new_pdb.append(new_line)

        elif line[0:6] == "HETATM":
            new_line = line[:60] + ("{bfac:6.2f}".format(bfac=default_val) ) + line[66:]
            new_pdb.append(new_line)

        else:
            new_pdb.append(line)

    new_pdb = "".join(new_pdb)
    with open(pdbOUT, 'w') as output_file:
        output_file.write(new_pdb)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def project_alnprofile(inpdat, mapDF, PDB):
    """
    Project vector of aln positions onto select protein
    Returns: dictionary with values for aligned resdiues in target PDB
    inpdat: input vector/1darray with length of alignment
    mapDF: mapping between alignment and protein positions
    PDB: target protein/PDB
    """
    # parse vector into PDB specific dict:
    res_dict = {}
    for i in range(len(inpdat)):
        if i in list(mapDF['aln']):
            current_res = int(list(mapDF[mapDF['aln']==i][PDB])[0]) + 1
            current_val = inpdat[i]

            res_dict[current_res] = current_val

    return res_dict



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def profile_byss(inputDF, pdblist, prefixOUT):
    """
    Projects sequence profiles onto the alignment
    """

    inpIdx = inputDF.index
    inputMAT = np.array(inputDF)

    sselm = pd.read_csv("../data/processed/consensus_aln_ss.txt", header=0, index_col=None, sep='\t')
    ssvec = ['b1', 'a1', 'b2', 'b3', 'a2', 'b4', 'a3', 'b5', 'a4', 'b6', 'a5']
    
    result = np.zeros(( len(pdblist), len(ssvec) ))
    for ix, i in enumerate(ssvec):
        current_sselm = sselm[sselm['element']==i]
        current_start = np.min(current_sselm['aln'])
        current_end = np.max(current_sselm['aln']) + 1
        current_mean = np.nanmean( inputMAT[:,current_start:current_end], 1)
        result[:,ix] = current_mean

    resultDF = pd.DataFrame(data=result, columns=ssvec, index=inpIdx)
    resultDF.to_csv("../data/processed/"+prefixOUT+"_aligned_ss.txt", header=True, index=True, sep='\t')
 
    # positions that don't align
    sel_nonaligned = np.ones(( np.shape(inputMAT)[1] ), dtype=bool)
    for i in sselm['aln']:
        sel_nonaligned[i] = False

    f_not = np.nanmean(inputMAT[:,sel_nonaligned], 1)
    f_aln = np.nanmean(inputMAT[:,sel_nonaligned==False], 1)
    fDF = pd.DataFrame({'not': list(f_not), 'aln': list(f_aln)})	

    print(fDF)
    fDF.to_csv("../data/processed/"+prefixOUT+"_aligned_avg.txt", header=True, index=False, sep='\t')







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_rmsf_data(pdblist, mapDF):
    """
    load rmsf data
    - computes rmsf over trajectory
    - averages over three replicas
    - stores median profile in dictionary
    """

    resdct = {}
    for PDB in pdblist:      
        print(PDB)
        sel_alnres = np.array(mapDF[PDB], dtype=int)
        run1 = get_rmsf("../gromacs/"+PDB+"_run1.xtc", "../gromacs/"+PDB+"_run1.gro", sel_alnres)
        run2 = get_rmsf("../gromacs/"+PDB+"_run2.xtc", "../gromacs/"+PDB+"_run2.gro", sel_alnres)
        run3 = get_rmsf("../gromacs/"+PDB+"_run3.xtc", "../gromacs/"+PDB+"_run3.gro", sel_alnres)
        current_run = np.vstack((run1, run2, run3))
        current_mean = np.mean(current_run, 0)
        resdct[PDB] = current_mean

    return resdct


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_rmsf(trajectory, topology, alnres=None):
    """
    get RMSF from trajectory
    uses aligned core residues for superposiion
    """

    reference = md.load(topology)
    top = reference.topology
    #table, bonds = top.to_dataframe()
    sel_top = top.select('protein and name CA')
    if len(alnres) > 0:
        sel_sup = sel_top[alnres]
    else:
        sel_sup = np.copy(sel_top)
        #print(sel_top)
        #print(len(sel_top))
        #print(table.iloc[sel_top])

    traj = md.load_xtc(trajectory, stride=2, top=topology)			# stride=2 bc computing resources
    traj.superpose(reference, frame=0, atom_indices=sel_sup)

    rmsf = md.rmsf(traj, reference, atom_indices=sel_top) #*10 #Convert from nm to angstroms
    #rmsf = md.rmsf(traj, reference)*10 #Convert from nm to angstroms

    return rmsf




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_distances(trajectory, topology, cmat):
    """
    get distances from trajectory
    IF CHANGING STRIDE, WATCH OUT FOR OUTPUT DIMENSIONS IN OTHER FUNCTIONS
    """

    #cmat = np.array([[10, 20], [30, 40]])

    reference = md.load(topology)
    traj = md.load_xtc(trajectory, stride=2, top=topology)			# stride=2 bc computing resources

    d = md.compute_distances(traj, atom_pairs=cmat)
    print(np.shape(d))

    return d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def aln2relpos(dfmap, poslist, orf):
    """
    maps aln position to input gene position
    dfmap: mapDF of corresponding positions
    poslist: input list of positions to query
    orf: input orf name
    """

    result = []
    for i in poslist:
        if i in list(dfmap['aln']):
            relpos = list(dfmap[dfmap['aln']==i][orf])[0]
            result.append(relpos)
        else:
            print("ERROR: not aln position")

    return result





#--------------------------------------------------------------------
def bond_clusters(distmat):
    """
    hierarchical bottom-up clustering of bonds based on input distance matrix
    """

    def average_linkage(clustermatrix, dmatIN, cluster1, cluster2):
        """
        compute average linkage of two clusters
        dmatIN: distance matrix on structure, eucledian distance
        cluster1, cluster2: id's to retrieve vals from dataIN
        """    
        if cluster1 in clustermatrix and cluster2 in clustermatrix:
            current_cluster1 = np.where(clustermatrix == cluster1)[0]
            current_cluster2 = np.where(clustermatrix == cluster2)[0]

            distance = 0
            counter = 0
            for i in current_cluster1:
                for j in current_cluster2:
                    current_distance = np.array([dmatIN[i,j], dmatIN[j,i]])
                    if np.any( np.isnan(current_distance) == False ):
                        distance += np.nanmax(current_distance)  
                        counter += 1
            if counter > 0:
                distance /= counter         #np.prod((len(current_cluster1), len(current_cluster2)))
            else:
                distance = np.nan
        else:
            distance = np.nan

        return distance


    def update_all(clustermatrix, dmat, dmatinit, cluster1, cluster2, prunelist, dendro):
        """
        update index table, cluster IDs as matrix indices
        intab: index table of cluster assignments
        cluster1, cluster2: clusters that are being merged
        """

        prunelist.append(int(cluster1))
        prunelist.append(int(cluster2))

        idx_cluster1 = np.where( clustermatrix == int(cluster1) )[0]
        idx_cluster2 = np.where( clustermatrix == int(cluster2) )[0]
        new_clustID = int(np.nanmax(clustermatrix)) + 1
        clustmat[idx_cluster1] = new_clustID
        clustmat[idx_cluster2] = new_clustID

        # update dendrogram
        dist = np.around( dmat[int(cluster1), int(cluster2)], 4 )
        size = len(clustmat[clustmat == new_clustID])
        dendro.loc[len(dendro)] = (cluster1, cluster2, new_clustID, dist, size)

        # recompute distance matrix
        nc,nr = dmat.shape
        new_dmat = np.zeros(( new_clustID+1, new_clustID+1 )) * np.nan  # copy to one added dimension
        new_dmat[0:nc,0:nc] = np.copy(dmat)


        for i in range(np.shape(new_dmat)[1]):
            if i not in prunelist:
                new_dmat[i, new_clustID] = average_linkage(clustmat, dmatinit, new_clustID, i)      # second cluster to iter

        new_dmat[np.array(prunelist, dtype=int),:] = np.nan
        new_dmat[:,np.array(prunelist, dtype=int)] = np.nan
        new_dmat[new_clustID, new_clustID] = np.nan         # check that new diagonal is nan
  
        return clustermatrix, new_dmat, prunelist, dendro



    ## INITIALIZE CLUSTERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dmat = np.copy(distmat) 
    clustmat = np.arange(np.shape(dmat)[0]) 
    dendro = pd.DataFrame(columns=["cluster1", "cluster2", "clusterNew", "distance", "size"]) 


    ## CLUSTERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    prunelist = []
    counter = 0
    dmatinit = np.copy(dmat)
    while len(list(set(clustmat))) > 0 and np.sum(np.isnan(dmat)==False) > 0:     # and counter < 2:      # (FOR DEBUGGING)
        min_i, min_j = np.unravel_index(np.nanargmin(dmat, axis=None), dmat.shape) 
        current_dist = dmat[min_i, min_j]
        #print(min_i, min_j, current_dist)
        
        # update
        clustmat, dmat, prunelist, dendro = update_all(clustmat, dmat, dmatinit, min_i, min_j, prunelist, dendro)

        counter += 1


    return dendro




#--------------------------------------------------------------------
def dendro_children(dndro, xid):
    """
    get children of cluster from dendrogram
    """

    if xid in list(dndro['clusterNew']):
        c1 = int(dndro[dndro['clusterNew']==xid]['cluster1'])
        c2 = int(dndro[dndro['clusterNew']==xid]['cluster2'])
        children = [c1, c2]
    else:
        children = None

    return children



#--------------------------------------------------------------------
def all_children(dndro, xid):
    """
    get all recursive children from dendrogram
    """

    list_int = [int(xid)] 
    list_int2 = [int(xid)]      # duplicate list to keep track of intermediates
    list_final = []
    while len(list_int) > 0:
        current_elem = list_int[-1]
        list_int.pop(-1)       # remove last element
        cc = dendro_children(dndro, current_elem)    
        if cc != None:
            list_int.append(cc[0])
            list_int.append(cc[1])
            list_int2.append(cc[0])
            list_int2.append(cc[1])
        else:
            list_final.append(current_elem)

    return list_final, list_int2



#--------------------------------------------------------------------
def clusters_max_dist(dndro, distmax):
    """
    get list of clusters with members
    """


    # get number of clusters from top 
    dndro = dndro[dndro['distance'] < distmax]
    cord = list(np.array(dndro['clusterNew'], dtype=int))[::-1]  
    list_clust = []
    list_deja = []
    for i in cord:
        if i not in list_deja:
            current_res, current_int2 = all_children(dndro, i)
            list_clust.append(current_res)
            list_deja.append(i)            
            list_deja += current_res
            list_deja += current_int2
   

    return list_clust




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def correlated_motions():
    """
    correlated motions through correlation of contact cluster distances
    """
    
    css = pd.read_csv("../data/processed/cluster_ss.txt", header=0, index_col=False, sep='\t')
    supermat = pickle.load( open("cluster_distances.pickle", 'rb') )
    print(np.shape(supermat))

    N = len(list_clustID)
    res = np.zeros(( N, N, len(pdblist) )) 
    for px, pdb in enumerate(pdblist):
        current_mat = supermat[:,:,px]
        
        for i in range( N ) :
            for j in range( i+1, N ): 
                current_cor = np.corrcoef(current_mat[i,:], current_mat[j,:])
                if current_cor[0,1] > 0.5:
                    res[i,j,px] = res[j,i,px] = 1
                   

    res = np.sum(res, axis=2) 
    resDF = pd.DataFrame(data=res, columns=list_clustID)
    resDF.insert(loc=0, column='clustID', value=list_clustID)
    resDF.to_csv("tmp.cm", header=True, index=True, sep='\t')		# len(list_clustID)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pdb_surface(PDB_file, pdbID, prefixOUT):
    """
    projects input vector onto pdb surface and writes to b-factor column
    wrangle with a combination of BioPython.PDBParser and Mdtraj 
    """

   
    # load and retrieve aggregation data
    agg_scores = pickle.load( open("../data/processed/agg_scores.pkl", 'rb') )
    current_agg = agg_scores[pdbID]
    
    # params
    theta_asa = 0.7
    kappa = 1
  
    # extract CA indices with PDBParser
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdbID, PDB_file)
    idx_CA = []
    for ax, atom in enumerate( structure[0]['A'].get_atoms() ):
        current_name = atom.get_name() 
        if current_name == "CA":
        	idx_CA.append(ax)
 

    # use Mdtraj  md.shrake_rupley with consistency with other analyses
    pdb = md.load(PDB_file)
    pdb_asa = md.shrake_rupley(pdb, mode='residue')[0]
    pdb_xyz = np.array(pdb.xyz)[0][idx_CA,:]

    exposed = pdb_asa > theta_asa
    exposed_xyz = pdb_xyz[exposed,:]
    exposed_agg = current_agg[exposed]
    exposed_idx = np.where(exposed)[0]
  
    score_dict = {}
    frame_score = np.zeros(( len(pdb_asa) ))
    for i in range(np.shape(exposed_xyz)[0]):
        query_xyz = exposed_xyz[i,:]
        query_dist = np.sum(np.sqrt(np.square(exposed_xyz - query_xyz)), 1)
        dist_weights = np.exp(-query_dist/kappa)
        weighted_score = np.sum(exposed_agg * dist_weights)
        frame_score[exposed_idx[i]] = weighted_score
        score_dict[exposed_idx[i]] = weighted_score
   

    color_bfac(PDB_file, score_dict, "../data/processed/"+pdbID+"_"+prefixOUT+".pdb")




if __name__ == '__main__':


    pdblist = [pdb.split('/')[-1].split('.')[0] for pdb in glob.glob('../data/pdb/*.pdb')]
    mapDF = pd.read_csv("../data/processed/c3718_relpositions.txt", header=0, index_col=None, sep='\t')
    #print(mapDF)


    # Compute frustration profiles
    frustdict = load_frustration_data(pdblist)
    frustDF, frustMAT = project_asa("/home/sebastian/foldingnets/code/c3718.aln", frustdict, "../data/processed/frustration_alnmat.txt")
    frustDF.to_csv("../data/processed/frustration_profiles.txt", header=True, index=True, sep='\t', na_rep="NA")


    frust_sd = np.nanstd(frustMAT, 0)
    frust_mean = np.nanmean(frustMAT, 0)
    frust_profile = pd.DataFrame({'pos': list(np.arange(len(frust_mean), dtype=int)+1), 'mean': list(frust_mean), 'sd': list(frust_sd)})
    frust_profile.to_csv("../data/processed/frustration_msaprofile.txt", header=True, index=True, sep='\t')

    bfac_frust = project_alnprofile(frust_sd, mapDF, 'YNL093W')
    color_bfac('../data/pdb/YNL093W.pdb', bfac_frust, "../data/processed/YNL093W_frustration.pdb")
    profile_byss(frustDF, pdblist, 'frustration')

    ###spectrum b, white_red, minimum=0, maximum=0.2
    ###spectrum b, white_0x006600



    # Compute rmsf profiles
    rmsfdict = load_rmsf_data(pdblist, mapDF)
    rmsfDF, rmsfMAT = project_asa("/home/sebastian/foldingnets/code/c3718.aln", rmsfdict, "../data/processed/rmsf_alnmat.txt")
    rmsfDF.to_csv("../data/processed/rmsf_profiles.txt", header=True, index=True, sep='\t', na_rep="NA")

    rmsf_sd = np.nanstd(rmsfMAT, 0)
    rmsf_mean = np.nanmean(rmsfMAT, 0)
    rmsf_profile = pd.DataFrame({'pos': list(np.arange(len(rmsf_mean), dtype=int)+1), 'mean': list(rmsf_mean), 'sd': list(rmsf_sd)})
    rmsf_profile.to_csv("../data/processed/rmsf_msaprofile.txt", header=True, index=True, sep='\t')


    bfac_rmsf = project_alnprofile(rmsf_sd, mapDF, 'YNL093W')
    color_bfac('../data/pdb/YNL093W.pdb', bfac_rmsf, "../data/processed/YNL093W_rmsf.pdb")
    profile_byss(rmsfDF, pdblist, 'rmsf')





    # DISTANCES, AVERAGES, PROFILES
    bpindex = pd.read_csv("../data/processed/cluster_members.txt", header=0, index_col=None, sep='\t')
    list_clustID = np.array(bpindex['clustID'])
    list_clustID = list_clustID[np.isnan(list_clustID)==False]		# remove non-cluster nans
    list_clustID = sorted(list(set(list_clustID)))					# non-redundant iter list
    print(len(list_clustID))


    supermat = np.zeros(( len(list_clustID), 15003, len(pdblist) ))		# ATTN: hard-coded dimension for trajectory length, watch out if chaging simulations!
    stdmat = np.zeros(( len(list_clustID), len(pdblist) ))		
    for px, pdb in enumerate(pdblist):
        print(pdb)
        cmat = np.array(bpindex[['A1', 'B1']])
        cmat[:,0] = aln2relpos(mapDF, cmat[:,0], pdb)
        cmat[:,1] = aln2relpos(mapDF, cmat[:,1], pdb)

        d1 = get_distances("../gromacs/"+pdb+"_run1.xtc", "../gromacs/"+pdb+"_run1.gro", cmat)
        d2 = get_distances("../gromacs/"+pdb+"_run2.xtc", "../gromacs/"+pdb+"_run2.gro", cmat)
        d3 = get_distances("../gromacs/"+pdb+"_run3.xtc", "../gromacs/"+pdb+"_run3.gro", cmat)
        dist = np.concatenate([d1, d2, d3], 0)

        for i in list_clustID:
            current_cluster = bpindex[bpindex['clustID']==i]
            if len(current_cluster) > 0:
                current_idx = np.array(current_cluster['i'], dtype=int) 	# 'i' is the index column in DF
                current_avg = np.nanmean(dist[:,current_idx], axis=1)		# average distance per frame over cluster
                current_std = np.nanstd(current_avg)						# std of cluster average profile
                supermat[int(i),:,px] = np.around(current_avg, 4)
                stdmat[int(i), px] = np.around(current_std, 4)

    np.savetxt("../data/processed/cluster_std.txt", stdmat)

    #pickle.dump(supermat, open("../data/processed/cluster_distances.pickle", 'wb'))
    #supermat = pickle.load( open("cluster_distances.pickle", 'rb') )


    # compute correlated motions from MD simulation data
    correlated_motions()



    # compute protein surface scores
    pdb_surface("../data/pdb/YNL093W.pdb", "YNL093W", "agg")

    traj = "YNL093W_07"
    PDB_file = "../data/CA/YNL093W_CA.pdb"
    trajectory = md.load("../data/trajectory/100ns_traj_1Tf_"+traj+".dcd", top=PDB_file)
    frame = trajectory[5000]
    frame.save_pdb("../data/processed/YNL093W_unfolded.pdb")
    pdb_surface("../data/processed/YNL093W_unfolded.pdb", "YNL093W", "unfolded_agg")




    ### !! Figure 3G !! change file name

    list_selclusters = [22, 28, 46, 71]

    resDF = pd.DataFrame(columns=['pdb', 'cluster', 'mean', 'std'])
    resmat = np.zeros(( len(pdblist), len(list_clustID) ))
    for px, pdb in enumerate(pdblist):
        current_mat = supermat[:,:,px]

        for ix, i in enumerate(list_clustID):			# list_selclusters
            current_cluster = current_mat[int(i),:]
            #current_cluster = np.random.choice(current_cluster, 1000, replace=False)

            current_mean = np.mean(current_cluster)
            current_std = np.std(current_cluster)
            resmat[px,ix] = current_mean
            resDF.loc[len(resDF)] = (pdb, i, current_mean, current_std)

            #for j in range(len(current_cluster)):
            #    current_dist = np.around(current_cluster[j], 4)

            #    resDF.loc[len(resDF)] = (pdb, i, j+1, current_dist)

    resDF.to_csv("tmp.selclust", header=True, index=False, sep='\t')
    
    np.savetxt("tmp.resmat", resmat)



