system_folders = ['system-setup', 'equilibrium', 'seeded-md', 'seeded-md/steering']


def __parse_seeds(seeds):
    if seeds == 'all':
        seeds = [int(folder.split('_')[1]) for folder in os.listdir() if folder.startswith('snapshot_') and '.pdb' not in folder]
        seeds.sort()
    elif isinstance(seeds, (list, tuple)):
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