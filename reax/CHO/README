## Example of using reactive force field for relaxation

### Installation

1. Get the source file from lammps `https://github.com/lammps/lammps`
2. Go to the `lammps/src` directory, enable `reax` and `reaxc` by `make yes-reax` and `make yes-user-reaxc`
3. Compile `lammps/lib/reax`. Go to the directory, and `make -f Makefile.gfortran`.
4. Return to the `lammps/src` directory, compile lammps. `make serial`
5. You should able to see an executable `lmp_serial`. This is the lammps main programe. 

### Run the relaxation of CHO

1. Generate the structure file by `pymatgen`. 
```python
from pymatgen.io.lammps.data import LammpsData
lmp = LammpsData.from_structure(structure)
lmp.write_file('data.CHO')
```
2. Run the program by `lmp_serial -in in.CHO`

### Final words
This tutorial is based on the LAMMPS reax CHO example. Please cite the following paper if you will use it in paper. 

Disclaimer:  Using these force fields for systems they 
have not been explicitly trained against may produce
unrealistic results.  Please see the README file in 
each subdirectory for more detailed information.

Hydrocarbon oxidation C/H/O:

     The follow information is reproduced from:

     "Chenoweth, K.; van Duin, A.C.T.; Goddard, W.A. 
     J. Phys. Chem. A 2008, 112, 1040-1053."

     - To obtain the H/C/O compound data required to 
     extend the hydrocarbon-training set, DFT 
     calculations were performed on the following systems: 
     (a) dissociation energies for various bonds 
     containing carbon, oxygen, and hydrogen.  The 
     ground state structure was obtained through 
     full geometry optimization.  Dissociation curves 
     were calculated by constraining only the bond length of 
     interest and re-optimization of the remaining 
     internal coordinates. Optimization was also performed 
     for the various angles and torsions associated 
     with C/H/O interactions.

