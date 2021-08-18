# -*- coding: utf-8 -*-
# pylint: disable-all
import numpy as np
from ase import units
import ase
from pprint import pprint
from matplotlib import cm
from useful_functions import get_fit_from_points
import matplotlib.pyplot as plt
from ase.utils import pickleload
import os
from ase.geometry import geometry
from pprint import pprint


class ForceExtrapolation:
    """
    Extrapolation based on the forces at a finite difference electric fields

    ...

    Attributes
    ----------

    findiffdb : ase db
        Finite difference database
    fields_to_choose : list
        Fields over which the finite differencing is done
    db_params : dict
        Parameters which will be fed into the ase db select option

    vibresults : dict
        Results stored from finite difference calculations
    dmudR : dict
        Change in dipole moment with finite diff dispacement
    dFdG : dict
        Change in the Forces along with application of Electric field

    Methods
    -------

    get_dmudR
        Gets the change in the dipole moment along z directions ( periodic calcs )
    get_dFdG
        Gets the change in the forces along application of the electric field


    """

    def __init__(self, vibresults, fields_to_choose, atomsIS, atomsFS,\
            direction, displacement):

        self.fields_to_choose = fields_to_choose
        self.atomsIS = atomsIS
        self.atomsFS = atomsFS
        self.vibresults = vibresults

        self.dmudR = {}  # change of dipole moment accompanying disp of atoms
        self.dFdG = {}  # Change in the force with the electric field applied
        self.q = {}  # Get the charge based on dFdG
        self.modes = []  # Modes for doing the dot product
        self.direction = direction
        self.displacement = displacement

        # Default get the information from the forces and the electic field
        self._get_reaction_path()

    def get_dmudR(self):
        for indice in self.vibresults:
            try:
                mu_p = self.vibresults[indice]['p'][
                    self.displacement][0.0]['dipole']
                mu_m = self.vibresults[indice]['m'][
                    self.displacement][0.0]['dipole']
            except KeyError:
                continue

            delta_mu = mu_p - mu_m
            dmudR = delta_mu / 2 / self.displacement
            self.dmudR[indice] = dmudR

    def get_dFdG(self):

        for indice in self.vibresults:

            fields_forces = {}
            fields = []

            # try:
            for field in self.vibresults[indice][self.direction][
                    self.displacement]:
                fields_forces[field] = self.vibresults[indice][self.direction][
                    self.displacement][field]['forces']
                fields.append(field)
            # except KeyError:
            # continue

            if not self.fields_to_choose:
                f_pos = sorted([ff for ff in fields if ff > 0.0])

            else:
                # fields are provided
                f_pos = self.fields_to_choose

            fmax = np.max(f_pos)
            fmin = np.min(f_pos)
            deltaf = fmax - fmin

            dFdG = ( -1 * fields_forces[fmax] + 8 * fields_forces[fmin] \
                        - 8 * fields_forces[-1*fmin] + fields_forces[-1*fmax] ) / 6.0 / 2 / (deltaf)

            try:
                self.dFdG.append(dFdG)
            except (AttributeError, KeyError):
                self.dFdG = []
                self.dFdG.append(dFdG)

    def get_q(self):
        """ Get the q by taking the dot product between the dFdG and the
        reaction mode
        """

        dFdG = []
        j = 0
        for i in range(3 * len(self.vibresults)):
            if (i + 1) % 3 == 0:
                # a z-component
                try:
                    differential = self.dFdG[j]
                except IndexError:
                    print('Missing data!')
                    continue
                dFdG.append([0, 0, differential[-1]])
                j += 1
            else:
                dFdG.append([0, 0, 0])
        dFdG = np.array(dFdG)
        mu_axes = dFdG.T[-1]
        # now dot product with the different modes available
        for index, mode in enumerate(self.modes):
            try:
                q = np.dot(mu_axes, mode)
            except ValueError:
                continue
            self.q[index] = q

    def _get_reaction_path(self):
        """Gets the reaction path from the atoms object
        """
        ## check if the atoms are on the same side of the unit cell
        cell = self.atomsIS.get_cell()  # same cell used in IS and FS hopefully
        # get the vector respresenting the difference of the two
        vector_all = self.atomsIS.get_positions() - self.atomsFS.get_positions(
        )
        indices = []
        for index in self.vibresults:
            indices.append(index)
        indices = np.array(indices)
        vector_all = np.array(vector_all)
        vectors = vector_all[indices]
        min_vec = []
        for v in vectors:
            vmin, vlen = geometry.find_mic(v, cell, pbc=True)
            min_vec.append(vmin)
        ravel_vec = np.ravel(min_vec)
        self.modes.append(ravel_vec / np.linalg.norm(ravel_vec))


class EigenModesHessian:
    """

    Collects the Hessian and Eignemodes for a given
    transition state calculation

    :param directory: Directory where the pickle files are stored
    :type directory: str
    :atoms: Atoms object which corresponds to the transition state
    :type atoms: atoms object

    """
    def __init__(self, atoms, directory, indices, dx, prefix):
        self.atoms = atoms  # atoms object for TS
        self.directory = directory  # directory where the calculation was done
        self.indices = indices  # indices making up the vibrations modes
        self.dx = dx  # displacement that the atoms had
        self.prefix = prefix  # prefix assigned to vibration file

        self.H = np.empty(
            (3 * len(self.indices), 3 * len(self.indices)))  # Dynamical matrix

        self.modes = []  # eigenmodes
        self.frequencies = []  # eigenfrequencies
        self.frequencies_cm = []  # Real frequencies in cm-1
        self.unit_vectors = [
        ]  # Unit vectors to perturb atoms along normal mode

        # get the Hessian
        self.Hessian()
        self.eigenmodes()

    def Hessian(self):
        # Get the equilibrium forces
        forces0 = pickleload(
            open(os.path.join(self.directory, self.prefix + '.eq.pckl'), 'rb'))
        # Gets the Hessian from the pickle file
        r = 0
        for i, index in enumerate(self.indices):
            for axes in 'xyz':  # Moved along all direction
                pickle_filename = f'{self.prefix}.{index}{axes}'
                forces_p = pickleload(
                    open(
                        os.path.join(self.directory,
                                     pickle_filename + '+.pckl'), 'rb'))
                forces_n = pickleload(
                    open(
                        os.path.join(self.directory,
                                     pickle_filename + '-.pckl'), 'rb'))
                self.H[r] = (forces_n -
                             forces_p)[self.indices].ravel() / 4. / self.dx
                r += 1

    def eigenmodes(self):
        self.H += self.H.copy().T
        m = self.atoms.get_masses()[self.indices]
        self.im = np.repeat(m**-0.5, 3)
        omega2, modes = np.linalg.eigh(self.im[:, None] * self.H * self.im)
        self.modes = modes
        s = units._hbar * 1e10 / np.sqrt(units._e * units._amu)
        self.frequencies = s * omega2.astype(complex)**0.5
        self.frequencies_cm = np.real(0.01 * units._e / units._c /
                                      units._hplanck * self.frequencies)
