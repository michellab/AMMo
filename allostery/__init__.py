"""A library for assessing allosteric modulation of proteins via an sMD/MSM protocol.
"""
import os as _os
import warnings as _warnings
_warnings.simplefilter('ignore', category=UserWarning)
_warnings.simplefilter('ignore', category=DeprecationWarning)
_warnings.simplefilter('ignore', category=Warning)


__all__ = ['analysis',
           'equilibrium',
           'msm',
           'setup',
           'steering',
           'utils']

# check for BioSimSpace and pytraj
try:
    import BioSimSpace
    del BioSimSpace
except ModuleNotFoundError:
    raise ModuleNotFoundError('BioSimSpace required: www.biosimspace.org')

try:
    import pytraj
    del pytraj
except ModuleNotFoundError:
    raise ModuleNotFoundError('pytraj required: https://amber-md.github.io/pytraj/latest/index.html')

# check for AMBERHOME
if 'AMBERHOME' not in _os.environ:
    print('An installation of AMBER is required: https://ambermd.org/. Please install AMBER and set the AMBERHOME environment variable')

# find GROMACS in PATH
_gmx = False
for loc in _os.environ['PATH'].split(':'):
    if _os.path.exists(loc):
        if 'gmx' in _os.listdir(loc):
            _gmx = True
            break
if not _gmx:
   print('A GROMACS installation is required: https://www.gromacs.org/. Please install GROMACS and include it in your PATH')

from . import equilibrium
from . import setup
from . import steering
from . import utils