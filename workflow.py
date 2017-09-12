import molflow.definitions as mf


wf = mf.WorkflowDefinition("uv-vis-spectrum")
__workflow__ = wf

# Functions
setup_forcefield = mf.Function(sourcefile='./prep.py',
                               funcname='setup_forcefield')

equilibrate = mf.Function(sourcefile='./production.py',
                          funcname='equilibrate')

sample = mf.Function(sourcefile='./production.py',
                     funcname='sample')

postprocess = mf.Function(sourcefile='./production.py',
                          funcname='postprocess')

# Inputs
mol = wf.add_input('molecule', type='mdt', description='Input small molecule.')

temperature = wf.add_input('temperature',
                           description='Conformational sampling temperature',
                           default=298.0,
                           type='float')


# DAG
prepped_mol = setup_forcefield(mol)
equilibrated = equilibrate(prepped_mol)
samples = sample(equilibrated)
qm_calcs, spectrum = postprocess(mol, samples)


# Outputs
wf.set_output(spectrum, 'spectrum.txt', type='File')
