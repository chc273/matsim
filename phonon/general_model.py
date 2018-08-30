from phonopy import Phonopy
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure
import matplotlib.pyplot as plt
import logging

__author__ = "Chi Chen"
__date__ = "Aug 30, 2018"

class PhononFromModel(object):
    """
    This class uses frozen phonon method to calculate the phonon dispersion curves and thermodynamic
    properties of a structure.
    The steps are as follows:
        1. Displaced structures are obtained from pristine structure via phonopy.
        2. The forces of each displaced structures are calculated by model
        3. The forces are passed to phonopy object and the force constant are produced
        4. Phonon is calculated once the force constant matrix is available.

    Args:
        model (model object): a model object that has a method "calculate_forces" that takes in a list of
            structures and return a list of N*3 forces, where N are the number of atoms in each structures
        structure (pymatgen structure): a pristine pymatgen structure
        atom_disp (float): small displacements in Angstrom

    """

    def __init__(self, model=None, structure=None, atom_disp=0.015, **kwargs):
        self.model = model
        self.structure = structure
        self.atom_disp = atom_disp
        self.qpoints, self.vertices = get_qpoints_and_vertices(self.structure)
        is_plusminus = kwargs.get('is_plusminus', True)
        is_diagonal = kwargs.get('is_diagonal', True)
        is_trigonal = kwargs.get('is_diagonal', False)
        supercell_matrix = kwargs.get('is_diagonal', None)
        ph_structure = get_phonopy_structure(self.structure)
        if supercell_matrix is None:
            supercell_matrix = np.eye(3) * np.array((1, 1, 1))
        self.phonon = Phonopy(unitcell=ph_structure, supercell_matrix=supercell_matrix)
        self.phonon.generate_displacements(distance=self.atom_disp,
                                      is_plusminus=is_plusminus,
                                      is_diagonal=is_diagonal,
                                      is_trigonal=is_trigonal)

        disp_supercells = self.phonon.get_supercells_with_displacements()
        # Perfect supercell structure
        init_supercell = self.phonon.get_supercell()
        # Structure list to be returned
        self.structure_list = [get_pmg_structure(init_supercell)]

        for c in disp_supercells:
            if c is not None:
                self.structure_list.append(get_pmg_structure(c))

        forces = self.model.calculate_forces(self.structure_list)
        self.phonon.set_forces(forces[1:])
        self.phonon.produce_force_constants()
        logging.info("Force constant produced") 

    def get_thermo_properties(self, mesh=[8, 8, 8], t_step=10, t_max=3000, t_min=0):
        self.phonon.set_mesh(mesh=mesh)
        self.phonon.set_thermal_properties(t_step=t_step,
                                           t_max=t_max,
                                           t_min=t_min)
        plt = self.phonon.plot_thermal_properties()
        return plt

    def get_bs_plot(self, points_per_line=50, figsize=(6, 4), ylim=[-1, 5]):
        bands = []
        bands = append_bands(bands, self.qpoints, points_per_line)
        self.phonon.set_band_structure(bands)
        q_points, self.distances, self.frequencies, eigvecs = self.phonon.get_band_structure()

        special_points = self.phonon._band_structure._special_points
        distance = self.phonon._band_structure._distance

        plt.figure(figsize=figsize)
        for j, (d, f) in enumerate(zip(self.distances, self.frequencies)):
            for i, freqs in enumerate(f.T):
                if i == 0 and j == 0:
                    plt.plot(d, freqs, "b-", lw=1)
                else:
                    plt.plot(d, freqs, "b-", lw=1)
        for sp in special_points:
            plt.axvline(x=sp, linestyle=':', linewidth=1, color='k')
        plt.axhline(y=0, linestyle=':', linewidth=1, color='k')

        plt.xticks(special_points, ["$\mathrm{%s}$" % i for i in self.vertices])
        plt.xlabel("Wave vector")
        plt.ylabel("Frequency (THz)")
        plt.xlim(0, distance)
        plt.ylim(ylim)
        plt.yticks(list(range(ylim[0], ylim[-1]+1)))
        plt.tight_layout()
        return plt


def append_band(bands, q_start, q_end, n_points):
    band = []
    for i in range(n_points+1):
        band.append(np.array(q_start) +
                    (np.array(q_end) - np.array(q_start)) / n_points * i)
    bands.append(band)

def append_bands(bands, qpoints, n_points_per_line=50):
    for i in range(len(qpoints)-1):
        append_band(bands, qpoints[i], qpoints[i+1], n_points_per_line)
    return bands


def get_qpoints_and_vertices(pmg_structure, symprec=0.01, angle_tolerance=5.0):
    """
    Obtain a list of high-symmetry q-points and the associated labels for
    phonon specturm calculation.
    Args:
        pmg_structure (Structure): A pymatgen structure object.
        symprec (float): Tolerance for symmetry finding.
        angle_tolerance (float): Angle tolerance for symmetry finding.
    """

    hsk = HighSymmKpath(pmg_structure, symprec=symprec,
                        angle_tolerance=angle_tolerance)
    qpoints = []
    vertices = []

    for qp in hsk.kpath["path"]:
        for v in qp:
            qpoints.append(hsk.kpath["kpoints"][v])
            vertices.append(v)

    return np.array(qpoints), vertices


if __name__ == "__main__":
    import numpy as np
    from io import StringIO
    import os
    from pymatgen.io.lammps.data import LammpsData
    from pymatgen.core import Structure
    from monty.tempfile import ScratchDir
    import shutil

    def read_forces(filename):
        print('current dir ', os.getcwd(), os.listdir())
        with open(filename, 'r') as f:
            lines = f.readlines()[9:]
        data = np.loadtxt(StringIO(''.join(lines)))
        return data

    def run_file(file):
        os.system('lmp_serial -in '+ file)
    def lammps_forces(structures, inputfile):
        forces = []
        with ScratchDir("."):
            shutil.copy(inputfile, 'in.forces')
            for s in structures:
                lmp_data = LammpsData.from_structure(s)
                lmp_data.write_file('data.dump')  # the structure is written to data.dump
                run_file('in.forces')
                forces.append(read_forces('force.data'))
        return forces

    class ForceModel(object):
        def __init__(self, filename='../data/in.forces_snap'):
            self.input_file = filename

        def calculate_forces(self, structures):
            forces = lammps_forces(structures, self.input_file)
            return forces

    s = Structure.from_file('./data/POSCAR')
    fm = ForceModel()
    model = PhononFromModel(model=fm, structure=s)
    print('model initialized')
    plt = model.get_bs_plot(points_per_line=50)
    plt.show()


