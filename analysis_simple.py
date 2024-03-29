#!/usr/bin/env python
# Calculate generalized correlation in MD trajectory
#
# Author: Jiahao Fan
# Date: 23-9-2022

# Load the python package
import os
from dynetan.toolkit import *
from dynetan.viz import *
from dynetan.proctraj import *
from dynetan.gencor import *
from dynetan.contact import *
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
import argparse

def format_seconds(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)

def DNAcorr(input_dir, output_dir, trajectory, topology, segids, cutoff, kneighb, windows, ncores): 
    # Number of windows created from full simulation.
    numWinds = int(windows) 
    # Cutoff for contact map (In Angstroms)
    cutoffDist = int(cutoff)
    # K neighbor number
    kNeighb = int(kneighb)
    # Create the object that processes MD trajectories.
    dnap = DNAproc()

    getNodeFromSel
    # Path where input files will searched and results be written.

    # topology file
    top_file = os.path.join(input_dir, topology)

    # trajectory file
    traj_file = [os.path.join(input_dir, trajectory)]
    # dcdFiles = [os.path.join(workDir, "decarboxylase.1.short.dcd")]

    # ligandSegID = "OMP"

    # Segment IDs for regions that will be studied.
    print("input segID argument!!!", segids)
    segIDs = segids

    # Residue name for solvent molecule(s)
    h2oName = ["SOL"]

    # Sampled frames per window
    numSampledFrames = 100

    # Number of sampled frames for automatic selection of solvent and ions.
    numAutoFrames = numSampledFrames*numWinds

    # Network Analysis will make one node per protein residue (in the alpha carbon)
    # For other residues, the user must specify atom(s) that will represent a node.
    customResNodes = {}

    # We also need to know the heavy atoms that compose each node group.

    usrNodeGroups = {}

    usrNodeGroups["SOL"] = {}
    usrNodeGroups["SOL"]["OH2"] = set("OH2 H1 H2".split())

    #################################
    ### Extra configuration

    # Minimum contact persistance (In ratio of total trajectory frames)
    contactPersistence = 0.5

    #################################
    ### Load info to object

    dnap.setNumWinds(numWinds)
    dnap.setNumSampledFrames(numSampledFrames)
    dnap.setCutoffDist(cutoffDist)
    dnap.setContactPersistence(contactPersistence)
    dnap.seth2oName(h2oName)
    dnap.setSegIDs(segIDs)
    dnap.kNeighb = kNeighb

    dnap.setNodeGroups(customResNodes)
    dnap.setNodeGroups(usrNodeGroups)

    #load the trajectory
    dnap.loadSystem(top_file,traj_file)

    dnap.checkSystem()

    dnap.selectSystem(withSolvent=False)

    dnap.prepareNetwork()

    # If your system is too large, you can turn off the "in memory" option, at a cost for performance.
    dnap.alignTraj(inMemory=True)

    # To speed-up the contact matrix calculation, a larger stride can be selected, at a cost for precision.
    dnap.findContacts(stride=1, verbose=0)

    dnap.filterContacts(notSameRes=True, notConsecutiveRes=False, removeIsolatedNodes=True, verbose=0)

    # output_dir = input_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    print(dnap.contactMatAll[0, :, :])
    np.savetxt(output_dir+"/contact.txt", dnap.contactMatAll[0, :, :], fmt='%.2e')
    plt.figure(1)
    plt.imshow(dnap.contactMatAll[0, :, :])
    plt.colorbar()
    plt.savefig(output_dir+"/contact_cut%s.png"%cutoffDist, dpi = 500)

    # We can calculate generalized correlaions in parallel using Python's multiprocessing package.
    start = datetime.now()
    print("Started at: %s\n" % str(start))
    dnap.calcCor(ncores=int(ncores))
    end = datetime.now()
    time_taken = format_seconds((end - start).seconds)
    print("- Total time: %s\n" % str(time_taken))
    for wind_index in range(numWinds):
        with open(output_dir+'/gencorrelation_k%s_window%s.txt'%(kNeighb, wind_index), "w") as f:
            np.savetxt(f, dnap.corrMatAll[wind_index], fmt='%1.15f')
    
    print("corr matrix shape",np.shape(dnap.corrMatAll))
    return dnap.corrMatAll

def plot_heatmap(correlation, input_dir, output_dir, windows):
    # output_dir = input_dir+"/gencorr"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    fig, ax = plt.subplots()
    for i in range(int(windows)):
        plt.figure(3)
        plt.xlim(0,correlation[i].shape[0])
        plt.ylim(0,correlation[i].shape[0])
        plt.imshow(correlation[i])
        if i == 0:
            plt.colorbar()
        fig.set_size_inches(8, 6, forward=True)
        plt.savefig(output_dir+"/gencorr_window%s.png"%(i+1), dpi = 500)
        
if __name__ == "__main__":

    #parse cmd arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input directory", default="./")
    parser.add_argument("--output", help="Output directory", default="./")
    #custom arguments
    parser.add_argument("--trajectory", help="Trajectory file")
    parser.add_argument("--topology", help="Referencce PDB file (must contain the same number of atoms as the trajectory)")
    parser.add_argument("--segids", nargs='*', default = ["SYSTEM"], help="segment ID in pdb file")
    parser.add_argument("--cutoff", default = 100, help="cutoff for determining contact matrix")
    parser.add_argument("--k", default = 7, help="cutoff for determining contact matrix")
    parser.add_argument("--windows", default = 1, help="number of windows to calculate the corr matrix")
    parser.add_argument("--ncores", default = 2, help="number of cores for parallel computation")
    

    args = parser.parse_args()
    if os.path.exists(os.path.join(args.output, "gencorrelation_k%s_window0.txt"%(args.k))):
        print("importing gencorr matrix..\n")
        tmp_list = [0]*int(args.windows)
        for wind_index in range(int(args.windows)):
            tmp_list[wind_index] = np.loadtxt(os.path.join(args.output, "gencorrelation_k%s_window%s.txt"%(args.k, wind_index)))
        gencorr = np.array(tmp_list)
    else:
        print("calculating gencorr...")
        gencorr = DNAcorr(args.input, args.output, args.trajectory, args.topology, args.segids, args.cutoff, args.k, args.windows, args.ncores)
        
    plot_heatmap(gencorr, args.input, args.output, args.windows)
    
    