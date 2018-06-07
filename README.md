# mcts_md
Python implementation for PACS_MD and MCTS_MD.
It uses GROMACS command for all of the manipulation　about MD simulation.
The required version is GROMACS 2018.1.
## USAGE
### preparation
Before running MD, We have to do some preparations which are adding ion, energy minimization, nvt equilibration and npt equilibration.
We need the .mdp files to do each manipulation of such preparations. The .mdp files here are intended only for use with this chignolin folding
sample. You have to set the parameters for each of your tasks.  
```
initialize.py() 
```
### PACS MD
In the PACS MD, the number of cycles is 100 and the number of parallele cascades is 5 in default configuration.
This program makes output files including trajectory file (merged_pacs.trr), rmsd to the target structure (rmsd_pacs.xvg).
```
python pacs_md.py
```
### MCTS MD
```
python mcts_md.py 
```

