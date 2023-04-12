# Install
```
pip install dynetan==1.2.0
```

# Example
```
cd test
python ../analysis_simple.py --trajectory traj_aligned.xtc --topology final.pdb
```

# Usage
In your own usage, you should provide a trajectory (.xtc et al.) and a topology file (.pdb et al.) with the same number of atoms, typically we only want to keep the protein in the trajectory files, so the topology file should also include only the protein. There may be some unexpected error if you include solvent information in the trajectory file. 

# Result
After executing the code, you will get the contact matrix
```
contact.txt, contact_cut100.png
```
This is the contact matrix, meaning "which two pairs will be used for calculating correlation", the default setting for determining contact is 100 Angstroms, meaning two Alpha Carbon within 100 Angstroms is considered to be "in contact". This is generally enough, if your system is too large and you still want to include the atoms far apart to be in contact, you can set parameter by "--cutoff".

And the correlation matrix
```
gencorrelation_k7.txt ./gencorr/gencorr_window1.png
```
The theory behind the generalized matrix can be found in Phys. Rev. E 69, 066138. There is a hyperparameter "k" in the algorithm, which defaultly is set to 7. You can change this hyperparameter by "-k". Another parameter is the number of windows, which defaultly is set to 1, meaning we use the whole trajectory to calculater the correlation matrix, you can set it to be larger than 1, meaning that the trajectory is split into several parts to calculate the corr matrix separately.

Another parameter you can set is "--ncores", which is the number of cores for parallel computation, but you should test this parameter in practive since not always faster with more cores.