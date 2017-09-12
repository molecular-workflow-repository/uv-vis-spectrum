__DOCKER_IMAGE__ = 'docker.io/autodesk/moldesign:moldesign_complete_0.8.0rc4'

from itertools import chain
import numpy as np
import moldesign as mdt
from moldesign import units as u


def equilibrate(mol, temperature):
    mol.set_energy_model(mdt.models.OpenMMPotential)
    mol.set_integrator(mdt.integrators.OpenMMLangevin, temperature=temperature*u.kelvin,
                       frame_interval=250*u.fs,
                       timestep=0.5*u.fs, constrain_hbonds=False, remove_rotation=True,
                       remove_translation=True, constrain_water=False)

    mol.minimize()
    mol.run(2.5 * u.ps)
    return mol


def sample(mol):
    traj = mol.run(5.0 * u.ps)
    return traj


def postprocess(qmmol, mdtraj):
    qmmol.set_energy_model(mdt.models.CASSCF, active_electrons=6,
                         active_orbitals=6, state_average=6, basis='sto-3g')

    post_traj = mdt.Trajectory(qmmol)
    for frame in mdtraj:
        qmmol.positions = frame.positions
        qmmol.calculate()
        post_traj.new_frame()

    return post_traj


def make_spectrum(qmtraj):
    qmmol = qmtraj.molecule

    wavelengths_to_state = []
    oscillators_to_state = []

    for i in range(1, len(qmmol.properties.state_energies)):
        wavelengths_to_state.append(
                (qmtraj.state_energies[:, i]-qmtraj.potential_energy).to('nm', 'spectroscopy'))
        oscillators_to_state.append([o[0, i] for o in qmtraj.oscillator_strengths])

    all_wavelengths = u.array(list(chain(*wavelengths_to_state)))
    all_oscs = u.array(list(chain(*oscillators_to_state)))
    counts, bins = np.histogram(all_wavelengths, weights=all_oscs, bins=50)

    outs = []
    for c,b in zip(counts, bins):
        outs.append("%s   %s" % (c,b))

    return '\n'.join(outs)
