from os import environ, listdir, path
from argparse import ArgumentParser


def __change_interpreter(file, interpreter):
    with open(file, 'r') as fl:
        contents = fl.readlines()

    contents[0] = f'#!{interpreter}\n'

    with open(file, 'w') as fl:
        fl.writelines(contents)
    return None


def __check_for_python3():
    home = '/'.join(path.realpath(__file__).split('/')[:-1])
    python_path = ''

    for loc in environ['PATH'].split(':'):
        if 'python3' in listdir(loc):
            python_path = f'{loc}/python3'
    
    if python_path == '':
        raise ValueError('python3 not found in $PATH')
    elif python_path != '/usr/bin/python3':
        for file in listdir(f'{home}/bin'):
            if not file.startswith('_'):
                print()
                __change_interpreter(f'{home}/bin/{file}', python_path)
    
    return None


def __set_home():
    # find location
    home = '/'.join(path.realpath(__file__).split('/')[:-1])

    # prepare file contents
    source_file = f'{home}/allostery.sh'
    contents = [f'export ALLOSTERY_HOME={home}\n',
                 'export PATH="$ALLOSTERY_HOME/bin:$PATH"\n\n',
                 'if [ -z "$PYTHONPATH" ]; then\n',
                 '  export PYTHONPATH="$ALLOSTERY_HOME"\n',
                 'else\n',
                 '  export PYTHONPATH="$ALLOSTERY_HOME:$PYTHONPATH"\n',
                 'fi\n']

    # write new file
    with open(source_file, 'w') as file:
        file.writelines(contents)
    print(f'ALLOSTERYHOME set to {home}')


    return home


def __set_location(location, home):
    if not path.exists(location):
        raise ValueError(f'{location} does not exist. Please specify an existing path using --location')

    contents = [f'home: {home}\n',
                f'location: {location}\n', 
                f'project: \n']
    with open(f'{home}/config', 'w') as file:
        file.writelines(contents)
    print(f'Default project location set to {location}')

    return None


def setup(location):
    print(f'{"#"*17}\n### ALLOSTERY ###\n{"#"*17}')
    
    # check for python 3
    print('Looking for python 3...', end='')
    __check_for_python3()
    print('done.')

    # set allostery home
    home = __set_home()

    # set default project location
    __set_location(location, home)

    print('-'*30)
    print(f'Setup complete. Please add "source {home}/allostery.sh" to your .bashrc')
    print('-'*30)


def __main__():
    parser = ArgumentParser(description='set up allostery')
    parser.add_argument('--location', type=str, default=f'{environ["HOME"]}/Documents', help=f'Default location of allostery projects. Default: {environ["HOME"]}/Documents')
    args = parser.parse_args()

    setup(args.location)


if __name__ == '__main__':
    __main__()