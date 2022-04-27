import yaml
import os

# get allostery settings
with open(f'{os.environ["ALLOSTERY_HOME"]}/config', 'r') as file:
    allostery_settings = yaml.load(file, Loader=yaml.FullLoader)

# get current project settings
__config_file = 'allostery_settings["project_location"]}/{allostery_settings["project"]}/.defaults/config'
if not os.path.exists(__config_file):
    __config_file = f'{os.environ["ALLOSTERY_HOME"]}/data/project_default'  # change to default if no project settings

with open(__config_file, 'r') as file:
    settings = yaml.load(file, Loader=yaml.FullLoader)
