import yaml
import os as os
import numpy as np

# get allostery settings
with open(f'{os.environ["ALLOSTERY_HOME"]}/config', 'r') as file:
    _allostery = yaml.load(file, Loader=yaml.FullLoader)

# get current project settings
__config_file = f'{_allostery["location"]}/{_allostery["project"]}/.defaults/config'
if not os.path.exists(__config_file):
    __config_file = f'{os.environ["ALLOSTERY_HOME"]}/data/project_default'  # change to default if no project settings

with open(__config_file, 'r') as file:
    _project = yaml.load(file, Loader=yaml.FullLoader)


_system_folders = ['system-setup', 'equilibrium', 'seeded-md', 'seeded-md/steering']


def __parse_seeds(seeds):
    if isinstance(seeds, (list, tuple)):
        if not all(isinstance(idx, int) for idx in seeds):
            raise TypeError('Seed indices must all be of type int')
        else:
            return seeds
    elif isinstance(seeds, str):
        if '-' in seeds:
            idx_range = [int(idx) for idx in seeds.split('-')]
            return np.arange(idx_range[0], idx_range[1]+1).tolist()
        elif ',' in seeds:
            return [int(idx) for idx in seeds.split(',')]
        else:
            return [int(seeds)]
    elif isinstance(seeds, int):
        return [seeds]
    else:
        raise TypeError('Seeds must be of type str, int, list or tuple')