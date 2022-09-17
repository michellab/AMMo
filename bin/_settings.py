import yaml
import os as os

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