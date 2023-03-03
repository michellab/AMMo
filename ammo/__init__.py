"""A library for assessing allosteric modulation of proteins via an sMD/MSM protocol.
"""
import os as _os
import warnings as _warnings


__all__ = ['analysis',
           'equilibrium',
           'msm',
           'setup',
           'steering',
           'utils']

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
   