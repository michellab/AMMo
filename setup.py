from os import environ, listdir, path, makedirs
from argparse import ArgumentParser


def __change_interpreter(file, interpreter):
    with open(file, 'r') as fl:
        contents = fl.readlines()

    contents[0] = f'#!{interpreter}\n'

    with open(file, 'w') as fl:
        fl.writelines(contents)
    return None


def __check_for_python3(interpreter):
    # find AMMO_HOME
    home = '/'.join(path.realpath(__file__).split('/')[:-1])

    import sys
    interpreter_now = sys.executable

    # check that interpreter exists
    if not path.exists(interpreter):
        interpreter_select=input(f'Warning: {interpreter} does not exist.\nDo you want use current interpreter {interpreter_now} as default (y/n)?\n')
        while interpreter_select.lower() not in ['y', 'n']:
            interpreter_select=input('Please enter y/n:\n')
        if interpreter_select.lower() == 'y':
            interpreter = interpreter_now
        else:
            raise SystemExit(f'Please specify a valid python interpreter using --interpreter')
    
    
    if interpreter_now != interpreter:
        print(f'The defult interpreter to run AMMo commands is {interpreter}, but the current interpreter is {interpreter_now}')
        interpreter_select=input(f'Do you want to change the default interpreter to {interpreter_now} (y/n/quit)?\n')
        while interpreter_select.lower() not in ['y', 'n', 'quit']:
            interpreter_select=input('Please enter y/n/quit:\n')
        if interpreter_select.lower() == 'quit':
            raise SystemExit('Exiting setup...')
        elif interpreter_select.lower() == 'y':
            interpreter = interpreter_now
    
    # change in scripts if needed
    if interpreter != '/usr/bin/python3':
        print(f'Changing interpreter to {interpreter_now} in scripts...')
        for file in listdir(f'{home}/bin'):
            if not file.endswith('.py'):
                __change_interpreter(f'{home}/bin/{file}', interpreter)
    
    return None


def __set_home():
    # find location
    home = '/'.join(path.realpath(__file__).split('/')[:-1])

    # prepare file contents
    source_file = f'{home}/ammo.sh'
    contents = [f'export AMMO_HOME={home}\n',
                 'export PATH="$AMMO_HOME/bin:$PATH"\n\n',
                 'if [ -z "$PYTHONPATH" ]; then\n',
                 '  export PYTHONPATH="$AMMO_HOME"\n',
                 'else\n',
                 '  export PYTHONPATH="$AMMO_HOME:$PYTHONPATH"\n',
                 'fi\n']

    # write new file
    with open(source_file, 'w') as file:
        file.writelines(contents)
    print(f'AMMO_HOME set to {home}')


    return home


def __set_location(location, home):
    while not path.exists(location):
        
        create = input(f'{location} does not exist.\nDo you want to create it (y/n)?\n')
        while create.lower() not in ['y', 'n']:
            create = input('Please enter y/n: ')
        if create.lower() == 'y':
            makedirs(location, exist_ok=True)
        else:
            raise SystemExit(f'Please specify an existing path using --location')

    contents = [f'home: {home}\n',
                f'location: {location}\n', 
                f'project: \n']
    with open(f'{home}/config', 'w') as file:
        file.writelines(contents)
    print(f'Default project location set to {location}')

    return None


def setup(interpreter, location):
    print(f'{"#"*34}\n### Allostery in Markov Models ###\n{"#"*34}')
    print()
    
    # check for python 3
    print('Looking for python 3...', end='\n')
    __check_for_python3(interpreter)
    print('done.')

    # set allostery home
    home = __set_home()

    # set default project location
    __set_location(location, home)

    print('-'*30)
    print(f'Setup complete. Please add "source {home}/ammo.sh" to your .bashrc')
    print('-'*30)
    print('Commands can be run using "ammo", for all available commands run "ammo -h"')


def __main__():
    parser = ArgumentParser(description='set up AMMo (Allostery in Markov Models)')
    parser.add_argument('--interpreter', type=str, default='/usr/bin/python3', help='Python interpreter to run AMMo commands. Requirements: python3 and a yaml installation. Default: /usr/bin/python3')
    parser.add_argument('--location', type=str, default=f'{environ["HOME"]}/Documents/ammo_projects', help=f'Default location of allostery projects. Default: {environ["HOME"]}/Documents/ammo_projects')
    args = parser.parse_args()

    setup(args.interpreter, args.location)


if __name__ == '__main__':
    __main__()