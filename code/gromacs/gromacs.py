import os, sys
import pandas as pd
import numpy as np
import subprocess
import glob
import shutil as sh



#--------------------------------------------------------------------
def prepare_pdb(PDB, dist):
    """
    prepare gromacs input files from pdb
    use python+subprocess to automate
    """

    # extract coordinate entries from PDB only
    cmd_atom = "grep '^ATOM' ../data/pdb/" + PDB + ".pdb > " + PDB + ".pdb"
    output = subprocess.run(cmd_atom, shell=True)

    # generate gromacs inputs
    cmd_gro = "gmx pdb2gmx -f " + PDB + ".pdb -o " + PDB + ".gro -p " + PDB + "_topol.top -i " + PDB + "_posre.itp -ff amber99sb-ildn -water tip3p -ignh"
    output = subprocess.run(cmd_gro, shell=True)

    # add bounding box
    cmd_box = "gmx editconf -f " + PDB + ".gro -o " + PDB + "_box.gro -c -d "+str(dist)+" -bt cubic"
    output = subprocess.run(cmd_box, shell=True)

    # solvate 
    cmd_solv = "gmx solvate -cp " + PDB + "_box.gro -o " + PDB + "_solv.gro -p " + PDB + "_topol.top -cs spc216.gro"
    output = subprocess.run(cmd_solv, shell=True)



#--------------------------------------------------------------------
def check_charge(PDB):
    """
    check net-charge of system and neutralize through adding ions
    this is done manually as not worth to automate
    """

    # check system net charge
    cmd_charge = "gmx grompp -f mdp/em.mdp -c " + PDB + "_solv.gro -p " + PDB + "_topol.top -o " + PDB + "_ions.tpr -po " + PDB + "_mdout.mdp -maxwarn 1"
    output = subprocess.run(cmd_charge, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    #output.communicate()[0].decode("utf-8")

    NP, NN = [0, 0]
    result = output.stdout.decode('utf-8').split('\n')
    for line in result:
        if "System has non-zero" in line:
            print(line)
            current_charge = float(line.split(':')[-1])

            if current_charge > 0:
                NP, NN = [0, np.abs(current_charge)]
            else:
                NP, NN = [np.abs(current_charge), 0]


    return int(NP), int(NN)


#--------------------------------------------------------------------
def neutralize_system(PDB, NP, NN):
    """
    check net-charge of system and neutralize through adding ions
    this is done manually as not worth to automate
    NP: number of positive ions
    NN: number of negative ions
    """

    cmd_neut = "gmx genion -s " + PDB + "_ions.tpr -o " + PDB + "_solv_ions.gro -p " + PDB + "_topol.top -pname NA -nname CL -np " + str(NP) + " -nn " + str(NN)
    output = subprocess.run(cmd_neut, shell=True)



#--------------------------------------------------------------------
def run_equilibration(PDB):
    """
    run minimization, equilibration simulations
    requires all input files in place
    """

    file_posre = PDB + "_posre.itp"
    sh.copy(file_posre, "posre.itp")


    # steepest-descent energy minimization
    cmd_em1 = "gmx grompp -f mdp/em.mdp -c " + PDB + "_solv_ions.gro -p " + PDB + "_topol.top -o " + PDB + "_em.tpr -po " + PDB + "_em_mdout.mdp"
    cmd_em2 = "gmx mdrun -s " + PDB + "_em.tpr -o " + PDB + "_em.trr -x " + PDB + "_em.xtc -c " + PDB + "_em.gro -e " + PDB + "_em.edr -g " + PDB + "_em.log -cpo " + PDB + "_em.cpt"
    output = subprocess.run(cmd_em1, shell=True)
    output = subprocess.run(cmd_em2, shell=True)
    
    # NVT equilibration for 1ns
    cmd_nvt1 = "gmx grompp -f mdp/nvt.mdp -c " + PDB + "_em.gro -p " + PDB + "_topol.top -o " + PDB + "_nvt.tpr -po " + PDB + "_nvt_mdout.mdp" # -r " + PDB + "_em.gro"
    cmd_nvt2 = "gmx mdrun -s " + PDB + "_nvt.tpr -o " + PDB + "_nvt.trr -x " + PDB + "_nvt.xtc -c " + PDB + "_nvt.gro -e " + PDB + "_nvt.edr  -g " + PDB + "_nvt.log -cpo " + PDB + "_nvt.cpt"
    output = subprocess.run(cmd_nvt1, shell=True)
    output = subprocess.run(cmd_nvt2, shell=True)

    # NPT equilibration for 1ns
    cmd_npt1 = "gmx grompp -f mdp/npt.mdp -c " + PDB + "_nvt.gro -t " + PDB + "_nvt.cpt -p " + PDB + "_topol.top -o " + PDB + "_npt.tpr -po " + PDB + "_npt_mdout.mdp -maxwarn 2" # -r " + PDB + "_nvt.gro"
    cmd_npt2 = "gmx mdrun -s " + PDB + "_npt.tpr -o " + PDB + "_npt.trr -x " + PDB + "_npt.xtc -c " + PDB + "_npt.gro -e " + PDB + "_npt.edr  -g " + PDB + "_npt.log -cpo " + PDB + "_npt.cpt"
    output = subprocess.run(cmd_npt1, shell=True)
    output = subprocess.run(cmd_npt2, shell=True)

    print("[equilibrations done]")


#--------------------------------------------------------------------
def run_simulation(PDB, runID):
    """
    run production simulations
    requires all input files in place
    """

    cmd_md1 = "gmx grompp -f mdp/md.mdp -c " + PDB + "_npt.gro -t " + PDB + "_npt.cpt -p " + PDB + "_topol.top -o " + PDB + "_md.tpr -po " + PDB + "_md_mdout.mdp -maxwarn 1"
    cmd_md2 = "gmx mdrun -s " + PDB + "_md.tpr -o " + PDB + "_"+runID+".trr -x " + PDB + "_"+runID+".xtc -c " + PDB + "_"+runID+".gro -e " + PDB + "_"+runID+".edr  -g " + PDB + "_"+runID+".log -cpo " + PDB + "_"+runID+".cpt"
    output = subprocess.run(cmd_md1, shell=True)
    output = subprocess.run(cmd_md2, shell=True)




#--------------------------------------------------------------------
def initialize_gromacs(PDB, dist):
    """
    run production simulations
    requires all input files in place
    """

    print(PDB)
    prepare_pdb(PDB, dist)
    pos, neg = check_charge(PDB)
    print(pos, neg)
    neutralize_system(PDB, pos, neg)
    # this requires interactiv keyboard input: enter 13 for solvent! 





if __name__ == '__main__':

    pdblist = [pdb.split('/')[-1].split('.')[0] for pdb in glob.glob('../data/pdb/*.pdb')]
    pdblist = sorted(pdblist)

    


    # PREPARE GROMACS INPUT | requires interactive keyboard input and manual check
    for i in pdblist:
        initialize_gromacs(i, 1.0)

    # after inspection, solvation box sizes were increased:
    initialize_gromacs("YLR293C", 1.5) 
    initialize_gromacs("YKL154W", 2.2) 
    initialize_gromacs("YFL038C", 1.5) 
    initialize_gromacs("YOR094W", 2.2) 
    initialize_gromacs("YDL192W", 1.5) 
    initialize_gromacs("YBR164C", 2) 
    initialize_gromacs("YMR138W", 2) 


    # RUN SIMULATIONS | this runs for a bit (...)
    for i in pdblist:
            run_equilibration(i)
            run_simulation(i, "run1")
            run_simulation(i, "run2")
            run_simulation(i, "run3")
