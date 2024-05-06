import os, sys
from sys import stdout
import numpy as np
import pandas as pd
import subprocess
import multiprocessing as mp
import glob
import time
import xml.etree.ElementTree as ET
import networkx as nx
import itertools as it
import pickle

from Bio import AlignIO
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import sbmOpenMM
import mdtraj as md






#--------------------------------------------------------------------
def relpositions(alnIN):
    """
    Generates dataframe of relative protein and alignment positions
    alnIN: input multiple seq alignment in clustal format of 16 sequences
    """

    alignment = AlignIO.read(open(alnIN), "clustal")
    list_id = []

    for ix, record in enumerate(alignment):
        current_record = record.id.split('_')[0]
        list_id.append(current_record)
    list_id = list_id[0:16]

    alignment = np.array(alignment)[0:16,:]
    aligned_positions = np.sum(alignment == '-', 0) == 0 

    # mapping from aln to seq    
    mapDF = pd.DataFrame()

    for i in range(np.shape(alignment)[0]):
        current_id = list_id[i]
        current_aln = alignment[i, :]
        current_aln[aligned_positions] = "*"
        current_pos = current_aln[current_aln != "-"]
        current_ref = np.where(current_aln == "*")[0] 
      
        mapDF[current_id] = list( np.where(current_pos == "*")[0] )
    mapDF['aln'] = current_ref
    mapDF.to_csv("../data/processed/c3718_relpositions.txt", header=True, index=False, sep='\t')

    return mapDF



#--------------------------------------------------------------------
def prepare_SBMs(inputDIR):
    """
    Wrapper to prepare PDBs and generate SBM input files with SMOG2
    !! requires SMOG2 install in path !!
    - runs smog_adjustPDB to reformat PDB for SMOG2
    - runs SMOG2 to generate coarse-grained CA contact map
    - deletes other files (not needed) 
    """


    pdblist = [pdb.split('/')[-1].split('.')[0] for pdb in glob.glob(inputDIR + '*.pdb')]
   
    for current_pdb in pdblist:
    
        pdbIN = inputDIR + current_pdb + '.pdb'
        pdbADJ = '../data/adj/' + current_pdb + '_adj.pdb'
        pdbSMOG = '../data/CA/' + current_pdb + '.contacts'

        cmd_adjPDB = 'smog_adjustPDB -removeH -removewater -i ' + pdbIN + ' -o ' + pdbADJ
        output = subprocess.run(cmd_adjPDB, shell=True)

        cmd_smog2 = 'smog2 -CA -dname tmpsmog -i ' + pdbADJ + ' -s ' + pdbSMOG
        output = subprocess.run(cmd_smog2, shell=True)

        tmp = glob.glob('tmpsmog.*')
        for i in tmp:
            os.remove(i)



#--------------------------------------------------------------------
def folding_temperature(PDB, temperature_range, steps, prefixOUT):
    """
    determines folding temperature Tf of input CA-SBM
    input PDB: pdbIN
    input contact file: contactIN
    output prefix: preOUT
    in part adapted from sbmOpenMM tutorials
    """

    pdbIN = '../data/pdb/' + PDB + '.pdb'
    contactIN = '../data/CA/' + PDB + '.contacts'
    preOUT = PDB

    sbmCAModel = sbmOpenMM.models.getCAModel(pdbIN, contactIN)

    #Create the output folder if it does not exists
    folderEnergy = '../data/energy'
    if not os.path.exists(folderEnergy):
        os.mkdir(folderEnergy)

    for temperature in temperature_range:
        run_simulation(PDB, sbmCAModel, temperature, steps, prefixOUT)


    def readOpenMMReporterFile(reporter_file):
        #Open the reporter file
        with open(reporter_file, 'r') as ef:
            #Store the lines 
            lines = ef.readlines()
            #Define a dictionary to store the data
            data = {}
            #read the header and create for each term an entry for the dictionary initialised to a list
            for r in lines[0].split(','):
                data[r.replace('#','').replace('"','').strip()] = []
            #read each value in the file and store it in the dictionary's lists.
            for i,r in enumerate(data):
                for line in lines[1:]:
                    #Convert the read string into a float for easy processing of numerical data
                    data[r].append(float(line.strip().split(',')[i]))
                
            #Convert each list into a numpy array
            for entry in data:
                data[entry] = np.array(data[entry])
            
        #return the created dictionary
        return data

    #Create the output folder if it does not exists
    folderName = '../data/heatCapacityData'
    if not os.path.exists(folderName):
        os.mkdir(folderName)

    auxdata = pd.DataFrame(columns=['file', 'temperature'])

    #Iterate over the temeperature range of interest.
    for temperature in temperature_range:
    
        #Define the name of the energy file for each temperature
        energy_file = folderEnergy+'/'+PDB+'_energy_'+str(temperature)+'_'+prefixOUT+'.data'
    
        #Read the energy data from each energy file
        simulationData = readOpenMMReporterFile(energy_file)
    
        #For easy reading we store the potential energy numpy array into a variable
        V = simulationData['Potential Energy (kJ/mole)']
    
        #We define the path name of the outputfile
        fileName = folderName+'/'+PDB+'_'+str(temperature)+'_'+prefixOUT+'.data'
    
        #Save the potential energies into a file using numpy.savetxt method
        np.savetxt(fileName, V, newline="\n")

        auxdata.loc[len(auxdata)] = (fileName, temperature)

    # run external pyWHAM install
    pywham_file = "../data/pywham/pywham_"+preOUT+".out"
    parse_pywham_xml(auxdata, 0.1, '../data/pywham/pywham_'+prefixOUT+'_'+PDB+'.xml', pywham_file)
    run_pywham('../data/pywham/pywham_'+prefixOUT+'_'+PDB+'.xml')
    
    temperature = []
    heat_capacity = []
    
    #Read the PyWham heat capacity output file
    with open(pywham_file, 'r') as hcf:
      
        #Iterate over the lines and store the values
        for line in hcf:
            ls = line.strip().split()            #line splitted by columns
            temperature.append(float(ls[0]))     #column 1
            heat_capacity.append(float(ls[1]))   #column 2
            
    result = pd.DataFrame({'T': temperature, 'HC': heat_capacity})
    result.to_csv("../data/processed/result_pywham_"+preOUT+".txt", header=True, index=False, sep='\t')

    idx_max = np.argmax(np.array(heat_capacity))

    T_max = temperature[idx_max]
    HC_max = heat_capacity[idx_max]

    return T_max




#--------------------------------------------------------------------
def parse_pywham_xml(simDF, binint, xmlOUT, whamOUT):
    """
    Parses Pywham input XML file
    simDF: intput data frame with info on simulations
    binint: interval bin
    xmlOUT: output xml file
    whamOUT: output result file
    """
   
    #binint = '0.1'

    data = ET.Element('WhamSpec')
    data.tail = '\n'
    data.text = '\n'
    general = ET.SubElement(data, 'General')
    general.tail = '\n'
    general.text = '\n'

    coords = ET.SubElement(general, 'Coordinates')
    cordname = ET.SubElement(coords, 'Coordinate')
    cordname.set('name', 'V')
    coords.tail = '\n'

    DCR = ET.SubElement(general, 'DefaultCoordinateFileReader')
    DCR.set('returnsTime', 'false')
    returnlist = ET.SubElement(DCR, 'ReturnList')
    returnlist.set('name', 'V')
    DCR.tail = '\n'

    binnings = ET.SubElement(general, 'Binnings')
    binning = ET.SubElement(binnings, 'Binning')
    binning.set('name', 'V')
    interval = ET.SubElement(binning, "Interval")
    interval.text = str(binint)
    binnings.tail = '\n'

    trajs = ET.SubElement(data, 'Trajectories')
    trajs.tail = '\n'
    trajs.text = '\n'
    templist = []
    for i in range(len(simDF)):
        current_file = simDF.iloc[i]['file']
        current_T = simDF.iloc[i]['temperature']
        templist.append(str(current_T))
        trajectory = ET.SubElement(trajs, "Trajectory")
        trajectory.set("T", str(current_T))
        EnFunc = ET.SubElement(trajectory, "EnergyFunction")
        EnFunc.text = 'V'
        CoordFiles = ET.SubElement(trajectory, "CoordinateFiles")
        CoordFile = ET.SubElement(CoordFiles, "CoordinateFile")
        CoordFile.text = current_file
        trajectory.tail = '\n'
    trajs.text = '\t'

    jobs = ET.SubElement(data, 'Jobs')
    jobs.tail = '\n'
    jobs.text = '\n'
    HeatCap = ET.SubElement(jobs, 'HeatCapacity')
    HeatCap.set('outFile', whamOUT)
    HeatCap.text = '\n'
    EnFunc = ET.SubElement(HeatCap, 'EnergyFunction')
    EnFunc.text = 'V'
    EnFunc.tail = '\n'
    Temps = ET.SubElement(HeatCap, 'Temperatures')
    Temps.text = ",".join(templist)
    Temps.tail = '\n'
    HeatCap.tail = '\n'

    b_xml = ET.tostring(data)
    with open(xmlOUT, "wb") as f:
        f.write(b_xml)



#--------------------------------------------------------------------
def run_pywham(xmlIN):
    """
    Wrapper to run pywham
    !! requires local install/version of pywham-1.2 !!
    """

    cmd_pywham = './pywham-1.2/bin/wham.py ' + xmlIN
    output = subprocess.run(cmd_pywham, shell=True)



#--------------------------------------------------------------------
def run_simulation(PDB, SBM, temp, steps, preOUT):
    """
    Run SBM simulations
    in part adapted from sbmOpenMM tutorials
    """
    folderEnergy = '../data/energy'

    #Define the name of the energy file for each temperature
    energy_file = folderEnergy+'/'+PDB+'_energy_'+str(temp)+'_'+ preOUT + '.data'

    #Define the integrator and context for the simulation at the defined temperature
    integrator = LangevinIntegrator(temp, 1/picosecond, 0.5*femtoseconds)
    simulation = Simulation(SBM.topology, SBM.system, integrator)

    #Set the initial coordinates
    simulation.context.setPositions(SBM.positions)
    
    #Add a SBM reporter that writes energies every 1 picosecond = 2000 steps (at 0.5 fs timestep).
    simulation.reporters.append(sbmOpenMM.sbmReporter(energy_file, 2000, sbmObject=SBM,
                                                  step=True, potentialEnergy=True, temperature=True))
    #Run each simulation for 1 ns = 2 million steps.
    start_time = time.time()
    simulation.step(steps)
    print("--- Finished simlation at T=%s in %s seconds ---" % (temp, (time.time() - start_time)))



#--------------------------------------------------------------------
def pywham2df(PATH):
    """
    In case of restart, this function collects the 
    results files and parses them into a DF
    PATH = dir + file prefix
    """

    resDF = pd.DataFrame(columns=['pdb', 'Tf', 'HC'])

    resfiles = glob.glob(PATH+'*.txt')
    for i in resfiles:
        orf = i.split('_')[-1].split('.')[0]
        data = pd.read_csv(i, header=0, index_col=False, sep='\t')

        idx_max = np.argmax(np.array(data['HC']))

        T_max = data['T'][idx_max]
        HC_max = data['HC'][idx_max]

        resDF.loc[len(resDF)] = (orf, T_max, HC_max)

    return resDF



#--------------------------------------------------------------------
def run_sbm_simulation(PDB, Tfold, PREFIX="01ns_traj_Tf_", sim_steps=200000000, idx_start=0):
    """
    Run SBM simulation at input folding temperature
    """

    pdbIN = '../data/pdb/' + PDB + '.pdb'
    contactIN = '../data/CA/' + PDB + '.contacts'
    preOUT = PDB

    sbmCAModel = sbmOpenMM.models.getCAModel(pdbIN, contactIN)

  
    integrator = LangevinIntegrator(Tfold, 1/picosecond, 0.5*femtoseconds)
    simulation = Simulation(sbmCAModel.topology, sbmCAModel.system, integrator)
    simulation.context.setPositions(sbmCAModel.positions)
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open('../data/CA/'+preOUT+"_CA.pdb", 'w'))


    #Create the output folder if it does not exists
    folderEnergy = '../data/energy'
   
    temperature = Tfold

    #"""
    for replica in range(70):       # hard coded here how many replicas to run
    
        #output file names
        energy_file = '../data/energy/'+PREFIX+PDB+'_'+str(idx_start+replica+1).zfill(2)+'.data'
        trajectory_file = '../data/trajectory/'+PREFIX+PDB+'_'+str(idx_start+replica+1).zfill(2)+'.dcd'
    
        #simulation setup
        integrator = LangevinIntegrator(temperature, 1/picosecond, 0.5*femtoseconds)
        simulation = Simulation(sbmCAModel.topology, sbmCAModel.system, integrator)
    
        #Set the initial coordinates
        simulation.context.setPositions(sbmCAModel.positions)
    
        #Add a SBM reporter that writes energies every 1 picosecond = 2000 steps (at 0.5 fs timestep).
        simulation.reporters.append(sbmOpenMM.sbmReporter(energy_file, 2000, sbmObject=sbmCAModel,
                                                  step=True, potentialEnergy=True, temperature=True))
        #Add a DCD reporter that writes coordinates every 1 picosecond = 2000 steps (at 0.5 fs timestep).
        simulation.reporters.append(DCDReporter(trajectory_file, 2000))
    
        #1000 ns = 2000 million steps. # 2000000000
        start_time = time.time()
        simulation.step(sim_steps)
    
        print("--- Finished simulation for replica %s in %s seconds ---" % (replica+1, (time.time() - start_time)))
    #"""



#--------------------------------------------------------------------
def traj_rmsd(traj_file, PDB_file):
    """
    RMSD by frame from trajectory relative to reference PDB
    """

    trajectory = md.load(traj_file, top=PDB_file)
    pdb = md.load(PDB_file)
    rmsd = md.rmsd(trajectory, reference)*10 #Convert from nm to angstroms
    
    print(rmsd)
    np.savetxt("tmp.rmsd", rmsd)



#--------------------------------------------------------------------
def traj_Qf(traj_file, PDB_file, contact_file):
    """
    Computes reaction coordinate Qf from trajectory
    code in part adapted from sbmOpenMM tutorials
    """

    sbmCAModel = sbmOpenMM.models.getCAModel(PDB_file, contact_file)


    #Get the indexes of the native contacts as a list
    contacts = [(c[0].index, c[1].index) for c in sbmCAModel.contacts.keys()]

    #Load input.pdb as the reference structure
    reference = md.load(PDB_file)

    #Compute the referenece distances
    ref_distances = md.compute_distances(reference, contacts) #Note that mdtraj uses nanometers as distance units

    #Define a 20% error for the calculation of the native contacts.
    error = 1.20

    #Load input.pdb as the toplogy and traj.dcd as the trajectory file
    trajectory = md.load(traj_file, top=PDB_file)
    
    #Calculate the native contact distances
    sim_distances = md.compute_distances(trajectory, contacts)
    
    #Evaluate if contacts are formed or not in each simulation frame
    formed_native_contacts = np.where(sim_distances <= ref_distances*error, 1, 0) 
    
    #Calculate the number of formed contacts in each frame
    n_native_contacts = np.sum(formed_native_contacts, axis=1) 
    
    #Calculate the fraction of native contacts formed
    Qf = n_native_contacts/len(contacts)    

    return Qf






#--------------------------------------------------------------------
def add_contacts(tfDF):
    """
    Adds number of contacts to Tf data frame
    """

    list_contacts = []

    for PDB in list(tfDF['pdb']): 

        pdbIN = '../data/pdb/' + PDB + '.pdb'
        contactIN = '../data/CA/' + PDB + '.contacts'
        sbmCAModel = sbmOpenMM.models.getCAModel(pdbIN, contactIN)

        contacts = [(c[0].index, c[1].index) for c in sbmCAModel.contacts.keys()]
        list_contacts.append( len(contacts) )

    tfDF['contacts'] = list_contacts

    return tfDF




#--------------------------------------------------------------------
def traj_dists(traj_file, PDB_file, contact_file):
    """
    RMSD by frame from trajectory relative to reference PDB
    """

    trajectory = md.load(traj_file, top=PDB_file)
    pdb = md.load(PDB_file)
    contacts = np.loadtxt(contact_file)
    contacts = contacts[:,[1,3]]
    print(contacts)
    contacts -= 1 		# python indexing 
    contacts = np.array(contacts, dtype=int)
    print(contacts)
    print(np.shape(contacts))

    d = md.compute_distances(trajectory, atom_pairs=contacts)
    print(d)
    print(np.shape(d))
    print(np.sum(d[0,:]))
    print(np.max(d[0,:]))
    print(np.sum(d[0,:]), np.sum(d[540,:]))

    rmsd = md.rmsd(trajectory, pdb)*10 #Convert from nm to angstroms
    
    print(rmsd)
    #np.savetxt("tmp.rmsd", rmsd)




#--------------------------------------------------------------------
def check_unfold(traj_file, PDB_file, contact_file):
    """
    get reaction coordinate Qf from trajectory
    code adapted from sbmOpenMM tutorials

    Detection of breakpoins:
    1. identify unfolding window, Qf > 0.3 and Qf < 0.7 (+/- 200, first full unfolding event with padding)
    2. bond breakage: first to exceed threshold of 2.5 within unfolding window
    """

    sbmCAModel = sbmOpenMM.models.getCAModel(PDB_file, contact_file)
    contacts = [(c[0].index, c[1].index) for c in sbmCAModel.contacts.keys()]

    # identify first full unfolding event
    current_qf = traj_Qf(traj_file, PDB_file, contact_file)

    reference = md.load(PDB_file)
    ref_distances = md.compute_distances(reference, contacts) #Note that mdtraj uses nanometers as distance units

    trajectory = md.load(traj_file, top=PDB_file)
    sim_distances = md.compute_distances(trajectory, contacts)

    #np.savetxt("tmp.dist", sim_distances)
    #d2 = sim_distances[:,468]

    N = len(sbmCAModel.contacts.keys())
    window_rm = 20		# not centered, so doesn't have to be uneven
    #result_break = np.zeros(( N )) * np.nan
    resultDF = pd.DataFrame(columns=['A', 'B', 'init', 'bp'])

    med_distances = np.zeros(( np.shape(sim_distances) )) * np.nan
    std_distances = np.zeros(( np.shape(sim_distances) )) * np.nan  
    qf_smoothed = np.zeros(( len(current_qf)-window_rm )) 
    all_init = np.mean( sim_distances[0:window_rm,:], axis=0)
    for i in range(np.shape(sim_distances)[0] - window_rm):
        med_distances[i,:] = np.mean( sim_distances[i:(i+window_rm),:], axis=0 ) 
        std_distances[i,:] = np.std(  sim_distances[i:(i+window_rm),:], axis=0 )
        qf_smoothed[i] = np.mean( current_qf[i:(i+window_rm)] )


    if len(np.where(qf_smoothed < 0.3)[0]) > 0:            # first < 0.3, if not smoothed then all
        current_lo = np.where(qf_smoothed < 0.3)[0][0]
    else:
        current_lo = np.where(current_qf < 0.3)[0][0]

    qf_smoothed = qf_smoothed[0:current_lo]        # only consider frames before
    current_up = np.where(qf_smoothed > 0.8)[0][-1]  # max before last > 0.8
    if current_up > 100:
        premax = np.argmax(qf_smoothed[(current_up-100):current_up])
        current_up = current_up - 100 + premax
    else:
        current_up = 0

    list_A = []
    list_B = []
    bp_rm = np.ones(( np.shape(med_distances)[1] )) * 20000     # non-break, default instead of missing values (remember!)
    for j in range(np.shape(med_distances)[1]):
        list_A.append( int(contacts[j][0]) )
        list_B.append( int(contacts[j][1]) )

        current_spltpts = np.where(med_distances[:,j] > 2.5)[0]         # threshold of 2.5 for bondbreakage
        current_splits = np.split( med_distances[:,j], current_spltpts )   

        if len(current_spltpts) > 0:           
            current_bp = np.copy(current_spltpts)           # bp's within unfold window!
            current_bp = current_bp[current_bp < (current_lo + 300)]
            current_bp = current_bp[current_bp > current_up ]
            if len(current_bp) > 0:
                #bp_rm[j] = current_bp[-1]      # last 'break'
                bp_rm[j] = current_bp[0]        # first 'break'
            

    resultDF = pd.DataFrame({"A": list_A, "B": list_B, "init": list(all_init), "bp": list(bp_rm)})
    resultDF.to_csv("tmp.bp", header=True, index=False, sep='\t')
    
    return resultDF 



#--------------------------------------------------------------------
def consensus_contacts(pdblist, theta=0):
    """
    Get all contacts between aligned positions
    !! contact files have to be in the hard coded path !!
    """

    # all contacts onto aln + stats
    alignment = AlignIO.read(open("../data/processed/c3718.aln"), "clustal")
    mapDF = pd.read_csv("../data/processed/c3718_relpositions.txt", header=0, index_col=None, sep='\t')

    N = np.shape(alignment)[1]
    cmap = dict()

    for i in pdblist:
        current_contact = np.loadtxt("../data/CA/"+i+".contacts")
        current_contact = current_contact[:,[1,3]]			# omit chain ID
        for j in range( np.shape(current_contact)[0] ):
            current_A = int(current_contact[j,0]) - 1		# aln pos starts at 0
            current_B = int(current_contact[j,1]) - 1 		# aln pos starts at 0
        
            if current_A in list(mapDF[i]) and current_B in list(mapDF[i]):
                aln_A = list(mapDF[mapDF[i]==current_A]['aln'])[0]
                aln_B = list(mapDF[mapDF[i]==current_B]['aln'])[0]

                if int(aln_A) < int(aln_B):
                    index_AB = str(aln_A) + "_" + str(aln_B) 
                else:
                    index_AB = str(aln_B) + "_" + str(aln_A)

                if index_AB in cmap.keys():
                    cmap[index_AB] += 1
                else:
                    cmap[index_AB] = 1

    cmap2 = dict()
    for i in cmap.keys():
        current_count = cmap[i]
        if current_count > theta:
            cmap2[i] = current_count

    cmat = np.zeros(( len(cmap2), 2 ), dtype=int)
    for ix, i in enumerate(cmap2.keys()):
        current_contact = i.split("_")
        cmat[ix,0] = int(current_contact[0])
        cmat[ix,1] = int(current_contact[1])

    return cmap2, cmat



#--------------------------------------------------------------------
def initial_dist(PDB_file, contact_file):
    """
    reference distances of input contacts
    """

    sbmCAModel = sbmOpenMM.models.getCAModel(PDB_file, contact_file)
    contacts = [(c[0].index, c[1].index) for c in sbmCAModel.contacts.keys()]


    reference = md.load(PDB_file)
    ref_distances = md.compute_distances(reference, contacts) #Note that mdtraj uses nanometers as distance units

    return ref_distances
    


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_dssp(fname):
    """
    reads pre-computed DSSP file into dataframe
    """
    #: max accessible surface area (square angstroms) for amino acids, from
    #: `Tien et al (2013) <https://doi.org/10.1371/journal.pone.0080635>`_.
    #: `MAX_ASA_TIEN[a]` is the max surface area for amino acid `a`.
    MAX_ASA_TIEN = {'A':129.0, 'R':274.0, 'N':195.0, 'D':193.0, 'C':167.0,
                'E':223.0, 'Q':225.0, 'G':104.0, 'H':224.0, 'I':197.0,
                'L':201.0, 'K':236.0, 'M':224.0, 'F':240.0, 'P':159.0,
                'S':155.0, 'T':172.0, 'W':285.0, 'Y':263.0, 'V':174.0}

    dssp = pd.DataFrame(columns=['Residue', 'Chain', 'AA', 'ASA', 'RSA', 'SS', 'X', 'Y', 'Z'])
    counter = 0
    df_counter = 0
    with open(fname) as f:
    #    next(f)
        for line in f:
            counter += 1

            if counter > 28:                    # remove header of the DSSP output
                res     = int(line[6:10])
                chain   = str(line[11])
                AA      = str(line[13])
                ASA     = int(line[35:38])
                RSA     = np.around(ASA / float(MAX_ASA_TIEN[AA]), 2)
                SS      = str(line[16])

                if SS.strip() in ['G', 'H', 'I']:
                    SStype = 'helix'
                elif SS.strip() in ['B', 'E']:
                    SStype = 'strand'
                else:
                    SStype = 'loop'

                X       = float(line[115:122])
                Y       = float(line[123:130])
                Z       = float(line[131:138])

                dssp.loc[df_counter] = [res, chain, AA, ASA, RSA, SStype, X, Y, Z]
                df_counter += 1

    return dssp

#coord = dssp[dssp['Residue'] == 3 ][['X', 'Y', 'Z']]
#print( dssp['SS'] )


#--------------------------------------------------------------------
def assign_sselements(mapDF, pdblist):
    """
    - gets DSSP secondary structure (ss) for all aligned positions
    - derives consensus secondary structure
    - assigns secondary structure elements [sstopo] corresponding to protein family
    - isolated positions are assigned as connecting loops rather than consensus elements
    """

    sstopo = ['b1', 'a1', 'b2', 'b3', 'a2', 'b4', 'a3', 'b5', 'a4', 'b6', 'a5']

    # contacts to all positions involved 
    cmap, cmat = consensus_contacts(pdblist)
    cpos = []
    for i in cmap.keys():
        current_contact = i.split("_")
        cpos.append(int(current_contact[0]))
        cpos.append(int(current_contact[1]))
    cpos = sorted(set(list(cpos)))

    dsspDF = pd.DataFrame(columns=pdblist)
    for pdb in pdblist:

        # aln positions to actual positions
        list_pos = []
        for i in cpos:
            current_pos = list(mapDF[mapDF['aln']==i][pdb])[0]
            list_pos.append(current_pos + 1)	# python indexing to PDB res id!
     
        dssp = load_dssp('../data/dssp/'+pdb+'.dssp')

        list_dssp = []        
        for i in sorted(list_pos):
            current_ss = list(dssp[dssp['Residue']==i]['SS'])[0]
            #print(pdb, i, current_ss)
            list_dssp.append(current_ss)
        dsspDF[pdb] = list_dssp

    topoDF = pd.DataFrame(columns=['aln', 'ss'])
    for i in range(len(dsspDF)):
        current_dssp = list(dsspDF.iloc[i])
        h = np.sum( np.array(current_dssp) == 'helix' )/len(current_dssp)
        s = np.sum( np.array(current_dssp) == 'strand' )/len(current_dssp)
        l = np.sum( np.array(current_dssp) == 'loop' )/len(current_dssp)

        pct = np.array([h, s, l])
        ix = np.argmax(pct)
        current_max = pct[ix]
        if np.sum(pct == current_max) == 1:
            current_ss = ['helix', 'strand', 'loop'][ix]
        else:
            if current_max == pct[0]:
                current_ss = 'helix'
            elif current_max == pct[1]:
                current_ss = 'strand'
            else:
                current_ss = 'ERROR'
        
        c = cpos[i]
        topoDF.loc[len(topoDF)] = (c, current_ss)

    # there are 2 isolated consensus position assigned as 'helix' 
    # that are changed to 'loop' for the purpose of automatically
    # assigning secondary structure elements of len>1 position
    for i in range(1, len(topoDF)-1):
        data_last = topoDF.iloc[i-1]
        data_curr = topoDF.iloc[i]
        data_next = topoDF.iloc[i+1]
        if data_curr['ss'] != 'loop':
            if data_last['ss'] == 'loop' and data_next['ss'] == 'loop':
                topoDF.iloc[i] = (topoDF.iloc[i]['aln'], 'loop')
            elif (data_curr['aln'] - data_last['aln'] > 10) and (data_next['aln'] - data_curr['aln'] > 10):
                topoDF.iloc[i] = (topoDF.iloc[i]['aln'], 'loop')

    resultDF = pd.DataFrame(columns=['aln', 'ss', 'element'])
    current_assign = 'loop'
    topo_idx = 0
    for i in range( len(topoDF) ):
        if i < len(topoDF) -1: 
            current_data = topoDF.iloc[i]
            next_data = topoDF.iloc[i+1]
        
            if current_data['ss'] != 'loop' and topo_idx < len(sstopo):
                #current_pos = list(mapDF[mapDF['aln']==current_data['aln']]['YNL093W'])[0]
                #print(i, current_data['aln'], current_pos, current_data['ss'], sstopo[topo_idx])
                resultDF.loc[len(resultDF)] = (current_data['aln'], current_data['ss'], sstopo[topo_idx])

                if (next_data['ss'] != current_data['ss']):
                    topo_idx += 1
                elif next_data['aln'] - current_data['aln'] > 6: 		# some slack between aln positions
                    topo_idx += 1

            else: 
                resultDF.loc[len(resultDF)] = (current_data['aln'], current_data['ss'], 'none')
        
        else:
            current_data = topoDF.iloc[i]
            last_data = topoDF.iloc[i-1]
            if current_data['ss'] == last_data['ss']:
                resultDF.loc[len(resultDF)] = (current_data['aln'], current_data['ss'], sstopo[topo_idx])
            else:
                resultDF.loc[len(resultDF)] = (current_data['aln'], current_data['ss'], 'none')

    resultDF.to_csv("../data/processed/consensus_aln_ss.txt", header=True, index=False, sep='\t')

    return resultDF
    


#--------------------------------------------------------------------
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
def check_Qf30(TRAJ, PDB, CONTACT):
    """
    checks if protein unfoldings during simulation
    QF30: frames with less than 30% of native contacts
    wrapper function
    """
    Qf = traj_Qf(TRAJ, PDB, CONTACT)
    Qf30 = np.sum(Qf < 0.3)/len(Qf)

    return Qf30


#--------------------------------------------------------------------
def breakpoint(PDB, pdblist, mapDF):
    """
    generates breakpoint DF: 
    - checks if trajetory includes unfolding event
    - if so, extracts break point for set of consensus contacts
    """
    #PDB = "YDL192W"
    globlist = glob.glob("../data/trajectory/100ns_traj_1Tf_"+PDB+"_*.dcd")
    rr = np.zeros(( len(globlist) ))

    cmap, cmat = consensus_contacts(pdblist)
    resultDF = pd.DataFrame({'alnA': list(cmat[:,0]), 'alnB': list(cmat[:,1])})
    initDF = pd.DataFrame({'alnA': list(cmat[:,0]), 'alnB': list(cmat[:,1])})
    cmat[:,0] = aln2relpos(mapDF, cmat[:,0], PDB)
    cmat[:,1] = aln2relpos(mapDF, cmat[:,1], PDB)
    np.savetxt("tmp.cmat", cmat, fmt='%i')	 # need to write to disk for read-in func


    counter_qf30 = 0
    for jx, j in enumerate(globlist):
    #for jx, j in enumerate(["../data/trajectory/100ns_traj_1Tf_YDL192W_07.dcd" ]): # "YDL192W_49" "YDL192W_53" "YDL192W_70" "YDL192W_22" "YDL192W_43" "YDL192W_25"
        current_alias = j.split('/')[3].split('.')[0].split('_')[3:]
        current_alias = "_".join(current_alias)

        Qf30 = check_Qf30(j, "../data/CA/"+PDB+"_CA.pdb", "../data/CA/"+PDB+".contacts")
        rr[jx] = Qf30

        if Qf30 > 0.01:         # min 1% of 'unfolded' frames
            print("UNFOLDING", PDB)
            counter_qf30 += 1
            bb = check_unfold(j, "../data/CA/"+PDB+"_CA.pdb", "tmp.cmat")  

            if counter_qf30 == 1:
                resultDF['A'] = bb['A']
                resultDF['B'] = bb['B']
                resultDF[current_alias] = bb['bp']
                initDF['A'] = bb['A']
                initDF['B'] = bb['B']
                initDF[current_alias] = bb['init']
            elif counter_qf30 > 1:
                resultDF[current_alias] = bb['bp']
                initDF[current_alias] = bb['init']
            else:
                print("ERROR")


    print(PDB, np.sum(rr > 0), np.sum(rr > 0)/len(rr))
        
    if os.path.exists("tmp.cmat"):
        os.remove("tmp.cmat")

    resultDF.to_csv("../data/breakpoint/bp_"+PDB+".txt", header=True, index=False, sep='\t')
    initDF.to_csv("../data/breakpoint/init_"+PDB+".txt", header=True, index=False, sep='\t')

    return resultDF



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def average_linkage_parallel( packaged_input ):
    """
    wrapper to run with multiprocessing library
    """
    
    clustermatrix = packaged_input[0]
    dataIN = packaged_input[1]
    cluster1 = packaged_input[2]
    cluster2 = packaged_input[3]

  
    if cluster1 in clustermatrix and cluster2 in clustermatrix:
        current_cluster1 = np.where(clustermatrix == cluster1)[0]
        current_cluster2 = np.where(clustermatrix == cluster2)[0]

        distance = 0
        for i in current_cluster1:
            for j in current_cluster2:
                distance += np.sqrt( np.square(dataIN[i] - dataIN[j]) ) 
        distance /= np.prod((len(current_cluster1), len(current_cluster2)))
 
    else:
        distance = np.nan

    return distance



#--------------------------------------------------------------------
def unfold_clusters(bpmat):
    """
    hierarchical bottom-up clustering of unfolding breakpoints
    """

    def compute_dmat(data):
        """
        compute initial distance matrix
        """

        N = len(data)
        dmat = np.zeros(( N, N )) * np.nan

        for i in range(N):
            for j in range(i+1, N):
                current_distance = np.sqrt( np.square(data[i] - data[j]) )
                dmat[i,j] = current_distance

        return dmat



    def average_linkage(clustermatrix, dataIN, cluster1, cluster2):
        """
        compute average linkage of two clusters
        dataIN: vector of breakpoints, eucledian distance
        cluster1, cluster2: id's to retrieve vals from dataIN
        """
    
        if cluster1 in clustermatrix and cluster2 in clustermatrix:
            current_cluster1 = np.where(clustermatrix == cluster1)[0]
            current_cluster2 = np.where(clustermatrix == cluster2)[0]

            distance = 0
            for i in current_cluster1:
                for j in current_cluster2:
                    distance += np.sqrt( np.square(dataIN[i] - dataIN[j]) ) 
            distance /= np.prod((len(current_cluster1), len(current_cluster2)))
 
        else:
            distance = np.nan

        return distance


    def update_all(clustermatrix, dmat, cluster1, cluster2, prunelist, dendro):
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
                new_dmat[i, new_clustID] = average_linkage(clustmat, data, new_clustID, i) 		# second cluster to iter

        #pool = mp.Pool(processes=20)
        #pool_set = np.where( [cluster not in prunelist for cluster in np.arange(new_clustID+1)] )[0]
        #pool_inputs = [(clustmat, data, new_clustID, i) for i in pool_set]
        #new_dmat[pool_set, new_clustID] = pool.map(average_linkage_parallel, pool_inputs)
        #pool.close()

        new_dmat[np.array(prunelist, dtype=int),:] = np.nan
        new_dmat[:,np.array(prunelist, dtype=int)] = np.nan
        new_dmat[new_clustID, new_clustID] = np.nan         # check that new diagonal is nan
  
        return clustermatrix, new_dmat, prunelist, dendro


    def parse_clustDF(clustmat, data):
        """
        compiles clusterDF with info on current clusters
        """

        clustDF = pd.DataFrame(columns=['min', 'max', 'mean', 'size', 'N'])
        for i in set(clustmat):
            current_members = np.where(clustmat==i)[0]
            current_data = data[current_members]
            clustDF.loc[len(clustDF)] = (np.min(current_data), np.max(current_data), np.mean(current_data), np.max(current_data)-np.min(current_data), len(current_data) )
        
        clustDF = clustDF.sort_values(by=['mean'], ascending=True)
        
        dvec = np.zeros(( len(clustDF) )) * np.nan
        current_min = np.array(clustDF['min'])
        current_max = np.array(clustDF['max'])
        dvec[0:(len(current_max)-1)] = current_min[1:len(current_min)] - current_max[0:(len(current_max)-1)]
        
        #for i in range(len(clustDF)-1):
        #    dvec[i] = clustDF.iloc[i+1]['min'] - clustDF.iloc[i]['max']
        clustDF['dist'] = list(dvec)
        clustDF['index'] = list( np.arange(len(clustDF)) + 1 )

        return clustDF



    def renumber_clusters(clustmat, data, clustDF):
        """
        renumber clusters for easier postprocessing
        """

        for i in set(list(clustmat)):
            current_cluster = clustmat == i
            current_ex = data[current_cluster][0]
            current_sel = (clustDF['min'] <= current_ex) & (clustDF['max'] >= current_ex)
            current_df = clustDF.loc[current_sel]
            current_id = int(list(current_df['index'])[0])
            clustmat[current_cluster] = current_id

        return clustmat


    ## INITIALIZE CLUSTERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    data = bpmat[np.isnan(bpmat) == False]
    clustmat = np.arange(len(data)) 
    dmat = compute_dmat(data)
    dendro = pd.DataFrame(columns=["cluster1", "cluster2", "clusterNew", "distance", "size"]) 


    ## CLUSTERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    prunelist = []
    counter = 0
    while len(list(set(clustmat))) > 1:  # and counter < 2:      # (FOR DEBUGGING)
        min_i, min_j = np.unravel_index(np.nanargmin(dmat, axis=None), dmat.shape) 
        current_dist = dmat[min_i, min_j]

        # update
        clustmat, dmat, prunelist, dendro = update_all(clustmat, dmat, min_i, min_j, prunelist, dendro)
        if len(list(set(clustmat))) < 5:
            clustDF = parse_clustDF(clustmat, data)
            print(clustDF)
            #if np.nanmin(np.array(clustDF['dist'])) > 20:	# dist btw clusters greater, break
            #    print(clustDF)
            #    break

        counter += 1

    clustmat = renumber_clusters(clustmat, data, clustDF)
    
    res = np.zeros(( len(bpmat) )) * np.nan
    res[np.isnan(bpmat) == False] = clustmat

    return res, dendro
    




#--------------------------------------------------------------------
def total_bpnan(theta=100):
    """
    This function filters out contacts that do not unfold systematically
    initially non-breaking bonds had a breakpoint of np.nan and
    systematically non-breaking bonds were filtered out here
    Change: non-breaking set to 20000, so filtered out by bp outside frame (> 10000)
    """


    # filter out contacts that do not systematically unfold as not informative for direct comparison
    for ix, i in enumerate(pdblist):
        bp = pd.read_csv("../data/breakpoint/bp_"+i+".txt", header=0, index_col=False, sep='\t', na_values="NA")
        bp2 = np.array(bp)[:,4:]
        if ix == 0:
            current_sum = np.sum(bp2 > 10000, axis=1)   #np.sum(np.isnan(bp2), axis=1)
        else:
            current_sum += np.sum(bp2 > 10000, axis=1)   #np.sum(np.isnan(bp2), axis=1)
    np.savetxt("../data/processed/totalnan_bp.txt", current_sum, fmt='%i')

    sel_nan = current_sum < theta
    
    #x2DF = xDF.iloc[sel_nan]
    #print(x2DF)

    return sel_nan




#--------------------------------------------------------------------
def bpstats(pdblist):
    """
    get suppl stats on bonds and breakpoints
    i,j: bond indeces
    A1/2, B1/2: aln indices of bond residues
    mean: mean breakpoint
    median: median breakpoint
    q80: 80pct quantile breakpoint

    """

    #cmap, cmat = consensus_contacts(pdblist)
    sel_nan = total_bpnan(100)

    dd = {}
    for pdb in pdblist:
        BPdf = pd.read_csv("../data/breakpoint/bp_"+pdb+".txt", header=0, index_col=False, sep='\t')
        BPdf = BPdf.iloc[sel_nan]
        current_cols = list(BPdf.columns)[4:]

        N = len(BPdf)
        M = len(current_cols)
        dmat = np.zeros(( N, N, M )) * np.nan

        for kx, k in enumerate(current_cols):
            data = np.array(BPdf[k]) 
            for i in range(N):
                for j in range(i+1, N):
                    current_distance = np.sqrt( np.square(data[i] - data[j]) )
                    dmat[i,j, kx] = current_distance    
        dd[pdb] = dmat
    ddmat = np.concatenate( [dd[x] for x in dd.keys()], axis=2)

    N, M, K = np.shape(dd[pdblist[0]])
    result = []
    resultDF = pd.DataFrame(columns=['i', 'j', 'A1', 'B1', 'A2', 'B2', 'mean', 'median', 'q80', 'nan', 'dist'])
    indexDF = pd.DataFrame(columns=['i', 'A1', 'B1'])
    dmat = np.zeros(( N, N )) * np.nan
    for i in range(N):
        A1 = BPdf.iloc[i]['alnA']   # reuse last BPdf from above, check
        B1 = BPdf.iloc[i]['alnB']
        indexDF.loc[len(indexDF)] = (i, A1, B1)
        for j in range(i+1, N): 
            A2 = BPdf.iloc[j]['alnA']
            B2 = BPdf.iloc[j]['alnB']
            distAB = np.array([np.sqrt(np.square(A1-A2)), np.sqrt(np.square(B1-B2))])
            if np.all(distAB <= 15): 
                current_mean  = np.nanmean(ddmat[i, j, :])
                current_median = np.nanmedian(ddmat[i, j, :])
                #current_q80 = np.percentile(ddmat[i,j,np.isnan(ddmat[i,j,:])==False], 80)
                current_q80 = np.percentile(ddmat[i,j, ddmat[i,j,:] < 10000], 80)
                current_nan   = np.sum( ddmat[i,j,:] > 10000) #np.sum(np.isnan(ddmat[i, j, :]))
                current_total = np.shape(ddmat)[2]

                tmp = np.zeros(( len(pdblist) ))
                for kx, k in enumerate(pdblist):
                    reference = md.load("../data/CA/"+k+"_CA.pdb")

                    cmat = np.array([[A1, A2], [B1, B2]])       # HERE DISTANCE WITHIN CHAIN !!!
                    cmat[:,0] = aln2relpos(mapDF, cmat[:,0], k)
                    cmat[:,1] = aln2relpos(mapDF, cmat[:,1], k)
                    ref_distances = md.compute_distances(reference, cmat) #Note that mdtraj uses nanometers as distance units
                    tmp[kx] = np.sum(ref_distances)

                resultDF.loc[len(resultDF)] = (i, j, A1, B1, A2, B2, np.around(current_mean, 4), current_median, current_q80, current_nan, np.around(np.mean(tmp), 4))
                dmat[i, j] = np.around(np.mean(tmp), 4)


    resultDF = resultDF.sort_values(by=['median'], ascending=True)
    resultDF.to_csv("../data/processed/bp_stats.txt", header=True, index=False, sep='\t')

    np.savetxt("tmp.dmat", dmat)
    indexDF.to_csv("../data/processed/bp_clusters.txt", header=True, index=False, sep='\t')

    return dmat, indexDF




#--------------------------------------------------------------------
def clusters_by_level(dndro, breakp, level):
    """
    get list of clusters with members
    """

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


    def all_children(dndro, xid):
        """
        get all recursive children from dendrogram
        """

        list_int = [int(xid)] #[1100]
        list_final = []
        while len(list_int) > 0:
            current_elem = list_int[-1]
            list_int.pop(-1)       # remove last element
            cc = dendro_children(dndro, current_elem)    
            if cc != None:
                list_int.append(cc[0])
                list_int.append(cc[1])
            else:
                list_final.append(current_elem)

        return list_final


    # get number of clusters from top
    cord = list(np.array(dndro['clusterNew'], dtype=int))[::-1]
    last = cord[0]
    list_clust = []
    for i in range(1,int(level)):
        last2 =  dendro_children(dndro, cord[i-1])
        list_clust += last2
        if cord[i-1] in list_clust:
            list_clust.pop(list_clust.index(cord[i-1]))

    #print(list_clust)
    list_result = []
    list_mean = []
    breakp = np.array(breakp)
    #print(len(breakp))
    for i in list_clust:
        current_res = all_children(dndro, i)
        #print(current_res)
        #print(breakp[current_res])
        list_result.append(current_res)
        current_mean = np.nanmean(breakp[current_res])
        current_max = np.nanmax(breakp[current_res])
        current_min = np.nanmin(breakp[current_res])
        #print(current_min, current_mean, current_max)
        list_mean.append(current_mean)

    amean = np.array(list_mean)
    sorder = np.argsort(amean)

    list_result_ordered = []
    for i in sorder:
        current_mean = amean[i]
        if current_mean < 20000:        # exclude non-break cluster
            list_result_ordered.append(list_result[i])

    return list_result_ordered





#--------------------------------------------------------------------
def bond_clusters(distmat):
    """
    hierarchical bottom-up clustering of bonds based on structural distance
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
def clusters_max_dist(dndro, bpidx, distmax=1):
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
   
    if len(bpidx) > 0:
        clust_id = np.zeros(( len(bpidx) )) * np.nan
        counter_id = 0
        for i in list_clust:
            current_idx = np.array(i, dtype=int)
            clust_id[current_idx] = counter_id
            counter_id += 1
            #print(bpidx.loc[i])
        bpidx['clustID'] = list(clust_id)

    bpidx.to_csv("../data/processed/cluster_members.txt", header=True, index=False, sep='\t', na_rep='NA')

    return list_clust, bpidx




#--------------------------------------------------------------------
def ccss(mapDF, pdblist, idxDF):
    """
    clustered contacts secondary structure assignment
    MOVE PRINT FILE UP TO SEPRATE INSTANCE EARLIER
    """

   
    sselm = assign_sselements(mapDF, pdblist)

    resDF = pd.DataFrame(columns=['clustID', 'ss_A1', 'ss_B1'])
    for i in list(set(list(idxDF['clustID']))):
        current_idxdf = idxDF[idxDF['clustID']==i]
        #print(current_idxdf)
        list_ss1 = []
        list_ss2 = []
        
        for j in range(len(current_idxdf)):
            current_A1 = current_idxdf.iloc[j]['A1'].item()
            current_B1 = current_idxdf.iloc[j]['B1'].item()

            ss_c1 = list(sselm[sselm['aln'] == current_A1]['element'])[0]
            ss_c2 = list(sselm[sselm['aln'] == current_B1]['element'])[0]
            list_ss1.append(ss_c1)
            list_ss2.append(ss_c2)

        current_ss1 = sorted(list(set(list_ss1)))
        current_ss2 = sorted(list(set(list_ss2)))

        if len(current_ss1) > 0 and len(current_ss2) > 0:
            # some cumbersome formatting
            if 'none' in current_ss1:
                current_ss1.pop(current_ss1.index('none'))
                current_ss1 = "/".join(current_ss1) + "*"
            else:
                current_ss1 = "/".join(current_ss1)
            if current_ss1 == "*":
                current_ss1 = "none"

            if 'none' in current_ss2:
                current_ss2.pop(current_ss2.index('none'))
                current_ss2 = "/".join(current_ss2) + "*"
            else:
                current_ss2 = "/".join(current_ss2)
            if current_ss2 == "*":
                current_ss2 = "none"
        
            resDF.loc[len(resDF)] = (i, current_ss1, current_ss2)
    
    
    resDF = resDF.sort_values(by=['clustID'], ascending=True)
    resDF.to_csv("../data/processed/cluster_ss.txt", header=True, index=False, sep='\t')

    return resDF
   



#--------------------------------------------------------------------
def representative_contacts(mapDF, pdblist, idxDF):
    """
    get most representative/populated of each contact cluster
    for plotting
    Attention: python indexing vs. indexing from external programm! 
    Attention: very few bondclusters are in non-aligned regions! 
    All output directly written to files for plotting! 
    """

    # 1. write native contacts into lists for check
    dict_contacts = {}

    for pdb in pdblist:
        native_contacts = np.loadtxt("../data/CA/"+pdb+".contacts")
        tmplist = []
        for j in range(np.shape(native_contacts)[0]):
            tmplist.append(str(int(native_contacts[j,1]))+"_"+str(int(native_contacts[j,3])))
        dict_contacts[pdb] = tmplist

    # 2. for each cluster (with bonds in aligned regions)_ write most representative bond
    resultDF = pd.DataFrame(columns=['A1', 'B1'])
    for ix, i in enumerate(sorted(list(set(idxDF['clustID'])))):
        current_cluster = idxDF[idxDF['clustID']==i]

        if len(current_cluster) > 0:
            tmpres = np.zeros(( len(current_cluster) ))
            for j in range(len(current_cluster)):
                current_coord = np.array(current_cluster.iloc[j][['A1', 'B1']]) 

                for pdb in pdblist:
                    rel_coord = aln2relpos(mapDF, current_coord, pdb)
                    rel_coord = np.array(rel_coord) + 1 # PYTHYON INDEXING ADJUSTMENT! 
                    check_coord = "_".join(sorted(list(np.array(rel_coord, dtype=str))))
                    native_contacts = dict_contacts[pdb]
                    if check_coord in native_contacts:
                        tmpres[j] += 1
            #print(i, tmpres)
            if np.sum(tmpres) > 0:
                best_rep = np.argmax(tmpres)
                resultDF.loc[len(resultDF)] = (current_cluster.iloc[best_rep][['A1', 'B1']])
    
    result = np.array(resultDF, dtype=int)
    np.savetxt("../data/processed/network_repclusterp.txt", result, fmt='%i')
    #print(resultDF)


    # 3. only aligned contacts
    Tf = pd.read_csv("../data/processed/Tf.txt", header=0, index_col=None, sep='\t')
    list_ints = []
    list_alns = []
    for pdb in list(Tf['pdb']):
        current_natives = dict_contacts[pdb]
        n_ints = len(current_natives)
        list_ints.append(n_ints)

        cmap, cmat = consensus_contacts(pdblist)    # all "aligned bonds"
        #cmat = np.array(idxDF[['A1', 'B1']])       # "aligned bonds" after filtering and clustering
        cmat[:,0] = aln2relpos(mapDF, cmat[:,0], pdb)
        cmat[:,1] = aln2relpos(mapDF, cmat[:,1], pdb)
        cmat = np.array(cmat, dtype=int) + 1 #python indexing! 

        current_count = 0
        for i in range(np.shape(cmat)[0]):
            check_coord = "_".join(sorted(list(np.array(cmat[i,:], dtype=str))))
            if check_coord in current_natives:
                current_count += 1
        list_alns.append(current_count)

    Tf['N_int'] = list_ints
    Tf['N_aln'] = list_alns

    Tf.to_csv("../data/processed/Tf_contacts2.txt", header=True, index=False, sep='\t')


    # 4. contact stats by secondary structure elemements
    sselm = assign_sselements(mapDF, pdblist)
    ssvec = ['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'a1', 'a2', 'a3', 'a4', 'a5', 'none']
    resultDF = pd.DataFrame(columns=ssvec)
    for pdb in pdblist:
        current_contact = np.loadtxt("../data/CA/"+pdb+".contacts")
        current_contact = current_contact[:,[1,3]] - 1       # omit chain ID + python indexing
        tmpres = np.zeros(( len(ssvec) ), dtype=int)

        for i in range(np.shape(current_contact)[0]):
            current_A = current_contact[i,0]
            current_B = current_contact[i,1]
            if current_A in list(mapDF[pdb]) and current_B in list(mapDF[pdb]):
                aln_A = list(mapDF[mapDF[pdb]==current_A]['aln'])[0]
                aln_B = list(mapDF[mapDF[pdb]==current_B]['aln'])[0]
                ss_A = list(sselm[sselm['aln'] == aln_A]['element'])[0]
                ss_B = list(sselm[sselm['aln'] == aln_B]['element'])[0]
                tmpres[ssvec.index(ss_A)] += 1
                tmpres[ssvec.index(ss_B)] += 1
        resultDF.loc[len(resultDF)] = tmpres
    resultDF.insert(loc=0, column='pdb', value=pdblist )
    resultDF.to_csv("../data/processed/contacts_aligned_ss.txt", header=True, index=False, sep='\t')


    # 5. distances in clusters
    distDF = pd.DataFrame(columns=['clustID', 'pdb', 'dist_res', 'dist_aa', 'N'])
    for ix, i in enumerate(sorted(list(set(idxDF['clustID'])))):
        current_cluster = idxDF[idxDF['clustID']==i]

        if len(current_cluster) > 0:
            current_coord = np.array(current_cluster[['A1', 'B1']]) 
            rel_coord = np.copy(current_coord)
            #res_rd = np.zeros(( len(pdblist) ))     # 3d distances of bonds
            #res_ad = np.zeros(( len(pdblist) ))     # distances between anchor AAs
            for px, pdb in enumerate(pdblist):
                rel_coord[:,0] = aln2relpos(mapDF, current_coord[:,0], pdb)
                rel_coord[:,1] = aln2relpos(mapDF, current_coord[:,1], pdb)
                rel_coord = np.array(rel_coord, dtype=int) + 1 # PYTHYON INDEXING ADJUSTMENT! 
                np.savetxt("tmp.cmat", rel_coord, fmt='%i')   
                rd = initial_dist("../data/CA/"+pdb+"_CA.pdb", "tmp.cmat")
                os.remove("tmp.cmat")       # remove tmp file

                # deal with edge cases
                aa_a = np.unique(rel_coord[:,0])
                if len(aa_a) > 1:
                    mat_a = np.array(list(it.combinations(np.unique(aa_a), 2)))
                aa_b = np.unique(rel_coord[:,1])
                if len(aa_b) > 1:
                    mat_b = np.array(list(it.combinations(np.unique(aa_b), 2)))
                if len(aa_a) > 1 and len(aa_b) > 1:
                    aa_coord = np.concatenate((mat_a, mat_b))
                elif len(aa_a) > 1 and len(aa_b) <= 1:
                    aa_coord = np.copy(mat_a)
                elif len(aa_a) <= 1 and len(aa_b) > 1:
                    aa_coord = np.copy(mat_b)
                np.savetxt("tmp.cmat", aa_coord, fmt='%i')   
                ad = initial_dist("../data/CA/"+pdb+"_CA.pdb", "tmp.cmat")
                os.remove("tmp.cmat")       # remove tmp file
               
                #res_rd[px] = np.mean(rd)
                #res_ad[px] = np.mean(ad)

                distDF.loc[len(distDF)] = (int(i), pdb, np.mean(rd), np.mean(ad), int(len(current_cluster)) )
    distDF.to_csv("../data/processed/clusters_avgdists.txt", header=True, index=False, sep='\t')





#--------------------------------------------------------------------
def consolidate_bp(pdblist, clusterlist, bpidx):
    """
    breakpoint consensus by cluster
    """
    
    tmp_df = pd.DataFrame(columns=['a', 'b'])
    list_df = []
    for pdb in pdblist:
        bp = pd.read_csv("../data/breakpoint/bp_"+pdb+".txt", header=0, index_col=False, sep='\t', na_values="NA")

        list_clustID = []
        resmat = np.zeros(( len(clusterlist), len(bp.columns)-4 ), dtype=int)     # -4 for annotation columns
        for ix, i in enumerate(clusterlist):
            current_clust = list(np.array(i, dtype=int))
            current_bp = bpidx.loc[current_clust]
            current_clustID = list(set(current_bp['clustID']))
            if len(current_clustID) == 1:
                current_clustID = int(current_clustID[0])
                list_clustID.append(current_clustID)
            else:
                print("ERROR")
                break
            list_idx = []
            for j in range(len(current_bp)):            # rechecking indices, superflous
                current_A = current_bp.iloc[j]['A1']
                current_B = current_bp.iloc[j]['B1']
                bp2 = bp[bp['alnA']==current_A]
                bp2 = bp2[bp2['alnB']==current_B]
                list_idx.append(list(bp2.index)[0])
            current_data = np.array( bp.iloc[list_idx][bp.columns[4:]] )
            
            current_consensus = np.zeros(( np.shape(current_data)[1] ))
            for k in range(np.shape(current_data)[1]):
                current_bond = current_data[:,k]

                #print(current_bond)
                dmat = np.zeros(( len(current_bond), len(current_bond) )) * np.nan
                for x in range(len(current_bond)):
                    for y in range(x+1, len(current_bond)):
                        dmat[x,y] =   np.sqrt(np.square(current_bond[x] - current_bond[y]))
                dndr = bond_clusters(dmat)
                #print(dndr)
                max_distance = np.max([20, (dndr.iloc[0]['distance']+1) ])
                clustlist, aa2 = clusters_max_dist(dndr,  tmp_df, distmax=max_distance)
                #print(clustlist)
                best_clust = 0
                best_len = 0
                for l in range(len(clustlist)):
                    if len(clustlist[l]) > best_len:
                        best_clust = l
                        best_len = len(clustlist[l])
                current_consensus[k] = np.median(current_bond[clustlist[best_clust]])
                #print(current_consensus[k])
        
            #current_consensus = np.array(np.around(np.median(current_data, axis=0), 0), dtype=int)
            resmat[ix,:] = current_consensus
 
        currentDF = pd.DataFrame(data=resmat, columns=list(bp.columns[4:]), index=list_clustID)
        list_df.append(currentDF)

    newDF = pd.concat(list_df, axis=1)
    newDF.insert(loc=0, column='clustID', value=list(currentDF.index) )
    
    return newDF



#--------------------------------------------------------------------
def clustpos_protein(indexDF, PDB):
    """
    retrieve cluster positions for input PDB
    """

    clustdict = {}

    for i in set(list(indexDF['clustID'])):
        current_idxdf = indexDF[indexDF['clustID']==i][['A1', 'B1']]
        if len(current_idxdf) > 0:
            vA = np.array(current_idxdf['A1'])
            vB = np.array(current_idxdf['B1'])
            vA = aln2relpos(mapDF, vA, PDB)
            vB = aln2relpos(mapDF, vB, PDB)

            print(i, vA, vB)

            clustdict[int(i)] = np.array( list(set(list(vA) + list(vB))) )

    return clustdict



#--------------------------------------------------------------------
def count_df(inplist):
    """
    get counts of elements in input list as df
    """

    countDF = pd.DataFrame(columns=['element', 'count'])
    for j in set(inplist):
        current_count = np.sum( np.array(inplist) == j )
        countDF.loc[len(countDF)] = (j, current_count)
    countDF = countDF.sort_values(by='count', ascending=False)

    return countDF



#--------------------------------------------------------------------
def classify_bpinit(allDF, bondssDF, start=0, topx=5):
    """
    Classify trajectories based on characteristic first breakpoints
    The categories have been defined based on instpection of the data 
    """
 
    list_columns = list(allDF.columns[1:])
    topx = 5
  
    # these are manually curated afer analysis
    list_b2 = ['a2', 'a2_b2', 'b2', 'b1_b2', 'a2_b3']    
    list_a5 = ['a5', 'a5_b4', 'a1_a5', 'a5_b1', 'a5_b5', 'a5_b6']   
    list_a3 = ['a3_a4', 'a3_none', 'a3', 'a2_a3', 'a3/b4', 'a3/b4_none']
    list_a1 = ['a1_none', 'a1_b2', 'a1', 'a1_b3']
    list_b5 = ['a4', 'a4_b5', 'a4_b4', 'b5', 'b5_none']


    res_assign = []
    res_not = []
    resultDF = pd.DataFrame(columns=['traj', 'protein', 'assign', 'type', 'class'])
    for i in range(len(list_columns)):
        current_column = list_columns[i]
        current_protein = current_column.split('_')[0]
        current_data = np.array(allDF[list_columns[i]])
        current_sort = np.argsort(current_data)
        current_top = current_sort[start:(start+topx)]
        current_bondss = bondssDF.iloc[current_top]

        dict_multi = {}
        dict_bondss = {}
        best_multi, best_multi_count = "none", 0
        best_bondss, best_bondss_count = "none", 0 
        for j in range(len(current_bondss)):
            tmp_bondss = current_bondss.iloc[j]
            tmp_bond = [tmp_bondss['ss_A1'].replace("*", ""),  tmp_bondss['ss_B1'].replace("*", ""), ]
            tmp_multi = "_".join(sorted(list(tmp_bond)))
            tmp_A1 = tmp_bond[0]
            tmp_B1 = tmp_bond[1]

            if tmp_multi in list(dict_multi.keys()):
                dict_multi[tmp_multi] += 1
                if dict_multi[tmp_multi] > best_multi_count:
                    best_multi_count = dict_multi[tmp_multi]
                    best_multi = tmp_multi
                    break
            else:
                dict_multi[tmp_multi] = 1

            if tmp_A1 in list(dict_bondss.keys()) and tmp_A1 != "none" and tmp_A1 != "b3":      # no 'b3' because it's buried and what links to it should cound (a2 or a5)
                dict_bondss[tmp_A1] += 1
                if dict_bondss[tmp_A1] > best_bondss_count:
                    best_bondss_count = dict_bondss[tmp_A1]
                    if best_bondss_count >= 2:
                        best_bondss = tmp_A1
            else:
                dict_bondss[tmp_A1] = 1

            if tmp_B1 in list(dict_bondss.keys()) and tmp_B1 != "none" and tmp_B1 != "b3":
                dict_bondss[tmp_B1] += 1
                if dict_bondss[tmp_B1] > best_bondss_count:
                    best_bondss_count = dict_bondss[tmp_B1]
                    if best_bondss_count >= 2:
                        best_bondss = tmp_B1
            else:
                dict_bondss[tmp_B1] = 1

        first_bondss = current_bondss.iloc[0]
        first_bond = [first_bondss['ss_A1'].replace("*", ""),  first_bondss['ss_B1'].replace("*", ""), ]
        first_bond = "_".join(sorted(list(first_bond))) 
 
        if best_multi != "none":
            res_assign.append(best_multi)
            current_assign = best_multi
            current_type = "multi"
        elif best_bondss != "none":
            res_assign.append(best_bondss)
            current_assign = best_bondss
            current_type = "best"
        else:
            res_assign.append(first_bond)
            current_assign = first_bond
            current_type = "first"

        if current_assign in list_b2:
            current_class = "a2/b2"
        elif current_assign in list_a5:
            current_class = "a5"
        elif current_assign in list_a3:
            current_class = "a3"
        elif current_assign in list_a1:
            current_class = "a1"
        elif current_assign in list_b5:
            current_class = "a4/b5"
        else:
            current_class = "other"
            #print(current_bondss)

        resultDF.loc[len(resultDF)] = (current_column, current_protein, current_assign, current_type, current_class)

    # parse data for network3D sankey network
    sankeyDF = pd.DataFrame(columns=['source', 'target', 'node_value', 'nodeID'])
    nodeID = 0
    for i in list(set(resultDF['protein'])):
        current_prot = resultDF[resultDF['protein']==i]
        for j in list(set(resultDF['class'])):
            current_class = current_prot[current_prot['class']==j]
            current_score = len(current_class)
            if current_score > 0:
                sankeyDF.loc[len(sankeyDF)] = (i, j, current_score, nodeID)
                nodeID += 1

    # summary stats of assignDF
    count_ass = count_df(list(resultDF['class']))
    count_ass['pct'] = count_ass['count'] / np.sum(np.array(count_ass['count']))


    resultDF.to_csv("../data/processed/assign_bpinit.txt", header=True, index=False, sep='\t')
    sankeyDF.to_csv("../data/processed/assign_sankey.txt", header=True, index=False, sep='\t')
    count_ass.to_csv("../data/processed/assign_counts.txt", header=True, index=False, sep='\t')

    return resultDF



#--------------------------------------------------------------------
def get_unfoldmat(allDF, bondssDF):
    """
    Compute relative unfolding curves
    allDF: all unfolding breakpoints
    bondssDF: secondary structure assignments for breakage clusters
    output nxmxk matrix
        n: secondary structure element
        m: order
        k: trajectory
    """


    ssdict = {}

    for i in range(len(bondssDF)):
        current_clustID = bondssDF.iloc[i]['clustID']
        ss1 = bondssDF.iloc[i]['ss_A1'].replace("*", "").split("/")
        ss2 = bondssDF.iloc[i]['ss_B1'].replace("*", "").split("/")
        ss = ss1 + ss2
        for j in ss:
            if j in list(ssdict.keys()):
                ssdict[j].append(current_clustID)
            else:
                ssdict[j] = list([current_clustID])


    list_ss = sorted(list(ssdict.keys()))
    print(list_ss)
    ufmat = np.zeros(( len(list_ss), len(allDF), len(list(allDF.columns)[1:]) ))

    for tx, traj in enumerate( list(allDF.columns)[1:] ):  
    
        current_data = np.array(allDF[traj])
        current_sort = np.argsort(current_data)

        # init
        current_bp = current_sort[0]
        for jx, j in enumerate(list_ss):
            if current_bp in ssdict[j]:
                ufmat[jx, 0, tx] += np.around(1/len(ssdict[j]), 3)
                #print(j)

        for i in range(1, len(current_sort)):
            current_bp = current_sort[i]
            current_uf = ufmat[:,i-1, tx]
            for jx, j in enumerate(list_ss):
                if current_bp in ssdict[j]:
                    current_uf[jx] += np.around(1/len(ssdict[j]), 3)
            ufmat[:,i, tx] = current_uf

    #print(ufmat)
    #print(np.shape(ufmat))
    #np.savetxt("tmp.ufmat", np.around(ufmat[2,:,:], 2))


    ssunfDF = pd.DataFrame(columns=['pdb', 'ss', 'pos', 'val'])
    for px, pdb in enumerate(pdblist):      # for each protein
        list_idx = []
        for jx, j in enumerate(list(allDF.columns[1:])):
            if pdb in j:
                list_idx.append(jx)
        ufmat_protein = np.mean( ufmat[:,:,list_idx], axis=2)

        for ix, i in enumerate(list_ss):   
            ufmat_ss = ufmat_protein[ix,:]

            for j in range(np.shape(ufmat_ss)[0]):
                current_val = np.around(ufmat_ss[j], 3)
                current_pos = j + 1
                current_ss = i
                current_prot = pdb
                #print(current_prot, current_ss, current_pos, current_val)
                ssunfDF.loc[len(ssunfDF)] = (current_prot, current_ss, current_pos, current_val)

#            res = np.zeros(( np.shape(ufmat)[1], len(pdblist) ))

#            res[:,px] = ufmat_protein[ix,:]

#        df = pd.DataFrame(data=res, columns=pdblist)
 #       df.insert(loc=0, column='pos', value=list(np.arange(np.shape(ufmat)[1])+1))
        #df.to_csv("../data/processed/ssunfold_"+i+".txt", header=True, index=False, sep='\t')
    ssunfDF.to_csv("../data/processed/ssunfold_byprotein.txt", header=True, index=False, sep='\t')


    bpinit = pd.read_csv("../data/processed/assign_bpinit.txt", header=0, index_col=False, sep='\t')
    ssunfDF2 = pd.DataFrame(columns=['start_ss', 'ss', 'pos', 'val'])
    aucDF = pd.DataFrame(columns=['start_ss', 'ss', 'auc', 'auc_sd'])
    list_columns = list(allDF.columns[1:])
    for px, ss in enumerate(list(set(bpinit['class']))):      # for each start ss
        current_sss = bpinit[bpinit['class']==ss]
        list_traj = list(current_sss['traj'])
        list_idx = []
        for traj in list_traj:
            current_idx = list_columns.index(traj)
            list_idx.append(current_idx)
        ufmat_protein = np.mean( ufmat[:,:,list_idx], axis=2)

        auc_protein = np.sum( ufmat[:,:,list_idx], axis=1)
        auc_mean = np.mean(auc_protein, axis=1)
        auc_sd = np.std(auc_protein, axis=1)

        for ix, i in enumerate(list_ss):   
            ufmat_ss = ufmat_protein[ix,:]
            current_ss = i
            aucDF.loc[len(aucDF)] = (ss, current_ss, auc_mean[ix], auc_sd[ix])
            for j in range(np.shape(ufmat_ss)[0]):
                current_val = np.around(ufmat_ss[j], 3)
                current_pos = j + 1
                current_prot = pdb
                ssunfDF2.loc[len(ssunfDF2)] = (ss, current_ss, current_pos, current_val)

    ssunfDF2.to_csv("../data/processed/ssunfold_bystart.txt", header=True, index=False, sep='\t')
    #print(ssunfDF2)
    aucDF.to_csv("../data/processed/ssunfold_bystart_auc.txt", header=True, index=False, sep='\t')
    print(aucDF)

    return ufmat





#--------------------------------------------------------------------
def compute_trajprops(allDF, pdblist):
    """
    computes various properties along trajectories:
    - aggregation surfaces 
    - Qf
    - radius of gyration
    Output to file in dataframes
    """

    cmap, cmat = consensus_contacts(pdblist)
    list_columnID = list(allDF.columns)[1:]

    agg_scores = pickle.load( open("../data/processed/agg_scores.pkl", 'rb') )
    theta_asa = 0.7
    kappa = 1
    result_agg = np.zeros(( 10000, len(list_columnID) ))  # trajectory length hard coded (!!)
    result_gyr = np.zeros(( 10000, len(list_columnID) ))  # trajectory length hard coded (!!)
    result_qf = np.zeros(( 10000, len(list_columnID) ))  # trajectory length hard coded (!!)


    for ix, i in enumerate( list_columnID ):     # omit first ID column
        current_pdb = i.split("_")[0]
        current_agg = agg_scores[current_pdb]
 
        traj = md.load("../data/trajectory/100ns_traj_1Tf_"+i+".dcd", top="../data/CA/"+current_pdb+"_CA.pdb")
        current_asa = md.shrake_rupley(traj, mode='residue')
        current_xyz = np.array(traj.xyz)
        #antiso = md.relative_shape_antisotropy(traj)
        current_gyr = md.compute_rg(traj)

        cmat_pdb = np.copy(cmat)
        cmat_pdb[:,0] = aln2relpos(mapDF, cmat_pdb[:,0], current_pdb)
        cmat_pdb[:,1] = aln2relpos(mapDF, cmat_pdb[:,1], current_pdb)  
        np.savetxt("tmp.cmat", cmat_pdb, fmt='%i')   # need to write to disk for read-in func
        current_qf = traj_Qf("../data/trajectory/100ns_traj_1Tf_"+i+".dcd", "../data/CA/"+current_pdb+"_CA.pdb", "tmp.cmat")
        os.remove("tmp.cmat")       # remove tmp file

        aggsurf = np.zeros(( len(traj) ))
        for frame in range(len(traj)):
            frame_asa = current_asa[frame,:]
            frame_xyz = current_xyz[frame,:,:]

            exposed = frame_asa > theta_asa
            exposed_xyz = frame_xyz[exposed,:]
            exposed_agg = current_agg[exposed]

            frame_score = 0
            for j in range(np.shape(exposed_xyz)[0]):
                query_xyz = exposed_xyz[j,:]
                query_dist = np.sum(np.sqrt(np.square(exposed_xyz - query_xyz)), 1)
                dist_weights = np.exp(-query_dist/kappa)
                weighted_score = np.sum(exposed_agg * dist_weights)
                frame_score += weighted_score
            aggsurf[frame] = frame_score

        result_agg[:,ix] = np.around(aggsurf, 2)
        result_gyr[:,ix] = current_gyr
        result_qf[:,ix] = current_qf

    resultDF = pd.DataFrame(data=result_agg, columns=list_columnID)
    resultDF.to_csv("../data/processed/surfaceAP_trajectory.txt", header=True, index=False, sep='\t')
    
    rgDF = pd.DataFrame(data=result_gyr, columns=list_columnID)
    rgDF.to_csv("../data/processed/rg_trajectory.txt", header=True, index=False, sep='\t')
     
    qfDF = pd.DataFrame(data=result_qf, columns=list_columnID)
    qfDF.to_csv("../data/processed/qf_trajectory.txt", header=True, index=False, sep='\t')
    



#--------------------------------------------------------------------
def bpintermediates_bycluster(allDF, clustering=False):
    """
    - cluster dendrograms of breakpoints
    - max difference between clusters
    - clustering: only need to be run if dendrograms need updating/dont exist
    """

   # CLUSTERING: used for 'intermediates'

    if clustering:
        for i in list(allDF.columns)[1:]:       # first columns is annotation
            current_clusters, current_dendro = unfold_clusters(np.array(allDF[i]))
            current_dendro.to_csv("../data/breakpoint/dendro/dendro_"+i+".txt", header=True, index=False, sep='\t')

    #intermediate breakpoint
    cmap, cmat = consensus_contacts(pdblist)
    pauseDF = pd.DataFrame(columns=['traj', 'start', 'end', 'mindist', 'last_bp', 'refold'])
    for tx, traj in enumerate( list(allDF.columns)[1:] ):  
    #list_traj = ['YDL192W_49','YDL192W_53','YDL192W_07','YDL192W_70','YDL192W_43','YDL192W_25']  
    #for tx, traj in enumerate( list_traj ):  
        current_dendro = pd.read_csv("../data/breakpoint/dendro/dendro_"+traj+".txt", header=0, index_col=None, sep='\t')    
        current_bp = np.array(allDF[traj])
        #np.savetxt("tmp.bp", current_bp)

        current_pdb = traj.split("_")[0]
        cmat_pdb = np.copy(cmat)
        cmat_pdb[:,0] = aln2relpos(mapDF, cmat_pdb[:,0], current_pdb)
        cmat_pdb[:,1] = aln2relpos(mapDF, cmat_pdb[:,1], current_pdb)  
        np.savetxt("tmp.cmat", cmat_pdb, fmt='%i')   # need to write to disk for read-in func
        current_qf = traj_Qf("../data/trajectory/100ns_traj_1Tf_"+traj+".dcd", "../data/CA/"+current_pdb+"_CA.pdb", "tmp.cmat")
        os.remove("tmp.cmat")       # remove tmp file
        #np.savetxt("tmp.qf", current_qf)
        #np.savetxt("../data/processed/qf_"+traj+".txt", current_qf)


        qf_smoothed = np.zeros(( len(current_qf) - 20 )) 
        for i in range(len(current_qf) - 20):
            qf_smoothed[i] = np.mean( current_qf[i:(i+20)] )

        if len(np.where(qf_smoothed < 0.3)[0]) > 0:
            current_lo = np.where(qf_smoothed < 0.3)[0][0]
        else:
            current_lo = np.where(current_qf < 0.3)[0][0]
        
        qf_refold = qf_smoothed[current_lo:]           # check for refolding events
        current_refold = current_lo + np.where(qf_refold > 0.8)[0]
        if len(current_refold):
            current_rf = current_refold[0]
        else:
            current_rf = 0

        qf_smoothed = qf_smoothed[0:current_lo]        # only consider frames before
        current_up = np.where(qf_smoothed > 0.8)[0][-1]  # max before last > 0.8
        if current_up > 100:
            premax = np.argmax(qf_smoothed[(current_up-100):current_up])
            current_up = current_up - 100 + premax
        else:
            current_up = 0

        
        best_pause = 0
        last_contact = 'none'
        for i in range(len(current_dendro)):
            cluster1 = int(current_dendro.iloc[i]['cluster1'])
            cluster2 = int(current_dendro.iloc[i]['cluster2'])

            current_res1, current_int2 = all_children(current_dendro, cluster1)
            current_memb1 = current_bp[current_res1]
            current_memb1 = current_memb1[current_memb1 < current_lo + 200]
            current_memb1 = current_memb1[current_memb1 > current_up]

            current_res2, current_int2 = all_children(current_dendro, cluster2)
            current_memb2 = current_bp[current_res2]
            current_memb2 = current_memb2[current_memb2 < current_lo + 200]
            current_memb2 = current_memb2[current_memb2 > current_up]
            
            if len(current_memb1) > 10 and len(current_memb2) > 10:         # for 'intermediates', at least 50 bps (~10%) evens before and after, no 20000ers
                min_dist = np.max([ np.min(current_memb1) - np.max(current_memb2), np.min(current_memb2) - np.max(current_memb1)])
                if min_dist > best_pause:
                    last_contact = np.min([np.max(current_memb1), np.max(current_memb2)])
                    best_pause = np.copy(min_dist)
        pauseDF.loc[len(pauseDF)] = (traj, current_up, current_lo, best_pause, last_contact, current_rf) 

    return pauseDF





#--------------------------------------------------------------------
def unfold_exdata(trajID, mapDF, pdblist):
    """
    Extracts data from exemplary unfolding trajectories for plotting
    """

    PDB = trajID.split("_")[0]
    traj_file = "../data/trajectory/100ns_traj_1Tf_"+trajID+".dcd"
    PDB_file = "../data/CA/"+PDB+"_CA.pdb"


    cmap, cmat = consensus_contacts(pdblist)
    current_sum = np.loadtxt("../data/processed/totalnan_bp.txt")
    sel_nan = current_sum < 100
    cmat = cmat[sel_nan,:]

    cmat[:,0] = aln2relpos(mapDF, cmat[:,0], PDB)
    cmat[:,1] = aln2relpos(mapDF, cmat[:,1], PDB)
    np.savetxt("tmp.cmat", cmat, fmt='%i')   # need to write to disk for read-in func


    sbmCAModel = sbmOpenMM.models.getCAModel(PDB_file, "tmp.cmat")
    contacts = [(c[0].index, c[1].index) for c in sbmCAModel.contacts.keys()]

    # output data
    current_qf = traj_Qf(traj_file, PDB_file, "tmp.cmat")
    trajectory = md.load(traj_file, top=PDB_file)
    sim_distances = md.compute_distances(trajectory, contacts)

    np.savetxt(trajID+"_distances.txt", sim_distances)
    np.savetxt(trajID+"_qf.txt", current_qf)












if __name__ == '__main__':


    pdblist = [pdb.split('/')[-1].split('.')[0] for pdb in glob.glob('../data/pdb/*.pdb')]
    print(pdblist)

    prepare_SBMs('../data/pdb/')
    mapDF = relpositions("../data/processed/c3718.aln")


    ## determine folding temperature Tf
    for i in pdblist:
        T_range = range(100,201,5)
        T_window = 10

        Tf1 = folding_temperature(i, T_range, 2000000, 'rough') 	
        Tf2 = folding_temperature(i, range(int(Tf1)-T_window,int(Tf1)+T_window+1), 20000000, 'refine') 	# 20 000 000
        Tf3 = folding_temperature(i, np.arange(int(Tf2)-2,int(Tf2)+2+0.5, 0.5), 20000000, 'fine')
        TfDF.loc[len(TfDF)] = (i, Tf3)

    TfDF = pywham2df("../data/processed/result_pywham_")
    TfDF.to_csv("../data/processed/Tf.txt", header=True, index=False, sep='\t')
    #TfDF = pd.read_csv("../data/processed/Tf.txt", header=0, index_col=False, sep='\t')


    TfDF2 = add_contacts(TfDF)
    TfDF2.to_csv("../data/processed/Tf_contacts.txt", header=True, index=False, sep='\t')



    ## GENERATE SBM SIMULATION DATA
    for i in pdblist:
        current_Tf = float(list(TfDF[TfDF['pdb']==i]['Tf'])[0])
        print(i, current_Tf)                                                                            # only 10ns contratry to file name :/
        run_sbm_simulation(i, current_Tf, PREFIX="100ns_traj_1Tf_", sim_steps=20000000,  idx_start=0)	# add index_start in python indexing




    ## analyse simulation data
    # gets all breakpoint DFs 
    for i in pdblist:
        bpDF = breakpoint(i, pdblist, mapDF)

    # get stats on breakpoints and structural distances
    struct_dmat, bpindex = bpstats(pdblist)
    #bpindex = pd.read_csv("../data/processed/bp_clusters.txt", header=0, index_col=None, sep='\t')
    #struct_dmat = np.loadtxt("tmp.dmat")

    struct_dendro = bond_clusters(struct_dmat)      # cluster bonds on structure
    list_clust, bpindex = clusters_max_dist(struct_dendro, bpindex, distmax=1)  # extracts clusters given distance cutoff
    #print(bpindex)

    allDF = consolidate_bp(pdblist, list_clust, bpindex)
    allDF.to_csv("../data/processed/allDF.txt", header=True, index=False, sep='\t')
    #allDF = pd.read_csv("../data/processed/allDF.txt", header=0, index_col=False, sep='\t')

    bondssDF = ccss(mapDF, pdblist, bpindex) 
    #bondssDF = pd.read_csv("../data/processed/cluster_ss.txt", header=0, index_col=False, sep='\t')

    assignDF = classify_bpinit(allDF, bondssDF)
    assignDF.to_csv("../data/processed/assign_bpinit.txt", header=True, index=False, sep='\t')

    ## representative contacts per cluster
    representative_contacts(mapDF, pdblist, bpindex)

    unfold_exdata("YNL093W_07", mapDF, pdblist)            # only one ex in final code

    ## compute relative unfolding curves
    ufmat = get_unfoldmat(allDF, bondssDF)      # output to files

    ## compute trajectory properties
    compute_trajprops(allDF, pdblist)      # output to files

    ## compute intermediates
    imdf = bpintermediates_bycluster(allDF, clustering=False)     # run clustering to generate/update dendrogram files
    imdf.to_csv("../data/processed/unfolding_intermediates.txt", header=True, index=False, sep='\t')

