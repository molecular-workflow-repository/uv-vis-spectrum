__DOCKER_IMAGE__ = 'docker.io/autodesk/moldesign:mdt_ambertools-0.8.0rc4'

import moldesign as mdt


def setup_forcefield(mol):
    """
    Create force field parameters for the chosen ligand
    """
    ff = mdt.create_ff_parameters(mol, charges='am1-bcc', baseff='gaff2')
    ff.assign(mol)
    return ff
