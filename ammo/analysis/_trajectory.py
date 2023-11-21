import os
from subprocess import run, DEVNULL
from pandas import read_csv, concat
from numpy import array, hstack, where
from matplotlib.pyplot import subplots as _subplots
from seaborn import kdeplot as _kdeplot
from pytraj import load_topology
from ._divergence import calculate_divergence
from ._hbonds import hbonds
from ammo.utils._utils import __check_cpptraj as _check_cpptraj, __get_trajectory as _get_trajectory, __random_name as _random_name

def compare_trajectories(traj, ref_traj, top=None, ref_top=None, n_bins=100, div_type='KL', filter=(20,20), filter_type='top', plot='plots', colors=('forestgreen', 'indigo'), pymol='analysis', residues=None, workdir=None, overwrite=False):
    """Compare two trajectories using divergence of hydrogen bond distances and phi, psi and chi1 dihedral angles
    
    Parameters
    ----------
    traj : str, pytraj.Trajectory
        the trajectory to analyse. If str, "top" must also be set

    ref_traj : str, pytraj.Trajectory
        the reference trajectory. If str, "ref_top" must also be set

    top : str
        topology file for the trajectory to be analysed

    ref_top : str
        topology file for the reference trajectory

    n_bins : int
        number of bins when calculating divergence

    div_type : str
        divergence type; "KL" or "JS"

    filter : [int, int], [float, float]
        top number (if int) or fraction (if float) of divergence/hbond frequency difference values to filter, in order of (hbonds, dihedrals); or the divergence cutoff(s). If None, all will be plotted/visualized
    
    filter_type : str
        way of filtering divergence. Allowed values are "top" (for top n or fraction values) and "cutoff" (for minimum divergence value(s))

    plot : string
        base file name for saving the plots of filtered distributions in the working directory

    colors : [str]
        colors for the trajectory and reference distributions, if plotting

    pymol : str
       base file name for saving a pymol session in the working directry. If None, no session will be saved

    residues : [int]
        list of residue indices (from 1) to carry out the comparison for. Useful if using wet trajectories or trajectories containing membranes. If None, all residues will be used

    workdir : str
        working directory for input/output files. If None, a directory will be created in the folder the script is run from

    overwrite : bool
        whether to overwrite the files that may be found in the working directory

    Returns
    -------
    
    """
    print('-'*30 + '\n-   Comparing trajectories   -\n' + '-'*30)
    # check that AMBERHOME is set
    _check_cpptraj()

    # fix divergence filter
    if filter is None:
        filter = [None, None]
    elif isinstance(filter, int) or isinstance(filter, float):
        filter = [filter, filter]

    # get working directory
    if workdir is None:
        workdir = f'./analysis_{_random_name()}'
        os.mkdir(workdir)
    else:
        if not os.path.exists(workdir):
            os.mkdir(workdir)

    # get trajectory file
    traj, top, remove_1 = _get_trajectory(traj, top)
    # get reference files
    ref_traj, ref_top, remove_2 = _get_trajectory(ref_traj, ref_top)

    # compute hbonds for trajectory and reference
    print(f'Finding H-bonds...', end='')
    all_hbonds = __find_hbonds(traj, top, ref_traj, ref_top, residues, workdir, overwrite)
    print('done.')
    # compute hbond distances for trajectory and reference
    print('Computing H-bond distances...', end='')
    traj_distances = __compute_hbond_distances(traj, top, all_hbonds.index, f'{workdir}/trajectory', overwrite)
    ref_distances = __compute_hbond_distances(ref_traj, ref_top, all_hbonds.index, f'{workdir}/reference', overwrite)
    #hbond_divergences = array([calculate_divergence(traj_distances[mask], ref_distances[mask], n_bins, div_type) for mask in all_hbonds.index])
    print('done.')

    # filter hbond divergence
    traj_distances, ref_distances, hbond_differences = __filter(traj_distances, ref_distances, all_hbonds['difference'], filter[0], filter_type)
 
    # compute all phi, psi and chi1 dihedral angles for trajectory and reference
    print('Computing torsional angles...', end='')
    traj_dihedrals = __compute_dihedrals(traj, top, residues, f'{workdir}/trajectory', overwrite)
    ref_dihedrals = __compute_dihedrals(ref_traj, ref_top, residues, f'{workdir}/reference', overwrite)
    dihedral_divergences = array([calculate_divergence(traj_dihedrals[mask], ref_dihedrals[mask], n_bins, div_type) for mask in traj_dihedrals.columns if mask in ref_dihedrals.columns])
    print('done.')

    # filter dihedral divergence
    traj_dihedrals, ref_dihedrals, dihedral_divergences = __filter(traj_dihedrals, ref_dihedrals, dihedral_divergences, filter[1], filter_type)

    # plotting
    fig, ax = __plot_distributions(traj_distances, ref_distances, hbond_differences, traj_dihedrals, ref_dihedrals, dihedral_divergences, colors, f'{workdir}/{plot}.png')
    
    # generate pymol
    if pymol is not None:
        pdb = __save_frame(traj, top, workdir)
        __generate_pymol(traj_distances.columns, all_hbonds['difference'], traj_dihedrals.columns, pdb, workdir, pymol)

    print('-'*30)
    return fig, ax


def __find_hbonds(traj, top, ref_traj, ref_top, residues, workdir, overwrite):
    if os.path.exists(f'{workdir}/hbonds.out') and not overwrite:
        all_hbonds = read_csv(f'{workdir}/hbonds.out', index_col=0)
        return all_hbonds

    if residues is None:
        hbond_mask = '*'
    else:
        hbond_mask = ','.join([f':{i}' for i in residues])
    traj_hbonds = hbonds(traj, hbond_mask, topology=top)
    ref_hbonds = hbonds(ref_traj, hbond_mask, topology=ref_top)
    # combine hbonds
    all_hbonds = concat([traj_hbonds, ref_hbonds], axis=1)
    all_hbonds.fillna(0, inplace=True)
    all_hbonds.columns = ['trajectory', 'reference']
    all_hbonds['difference'] = all_hbonds['trajectory'] - all_hbonds['reference']
    # save hbonds
    all_hbonds.to_csv(f'{workdir}/hbonds.out')

    return all_hbonds


def __compute_hbond_distances(traj, top, masks, fname, overwrite):
    input_file = f'{fname}_distances.in'
    output_file = f'{fname}_distances.out'

    if os.path.exists(output_file) and not overwrite:
        df = read_csv(output_file, delim_whitespace=True).iloc[:,1:]
        return df

    #parsed_masks = [' '.join([':'+atom.split('_')[1] for atom in mask.split('-')[:2]]) for mask in masks]
    parsed_masks = []
    for mask in masks:
        atom_masks = mask.split('-')
        residues = [mask.split('_')[1].split('@')[0] for mask in atom_masks[:2]]
        atom1 = atom_masks[0].split('_')[1].split('@')[1]
        atom2 = atom_masks[2]
        parsed_masks.append(f':{residues[0]}@{atom1} :{residues[1]}@{atom2}')
    with open(input_file, 'w') as file:
        file.writelines([f'parm {top}\n',
                         f'trajin {traj}\n'])
        for title, mask in zip(masks, parsed_masks):
            file.writelines([f'distance {title} {mask} out {output_file}\n'])
        file.writelines(['go\nquit\n'])
    
    run(['cpptraj', '-i', input_file], stdout=DEVNULL)

    df = read_csv(output_file, delim_whitespace=True).iloc[:,1:]
    return df


def __compute_dihedrals(traj, top, residues, fname, overwrite):
    input_file = f'{fname}_dihedrals.in'
    output_file = f'{fname}_dihedrals.out'

    if os.path.exists(output_file) and not overwrite:
        df = read_csv(output_file, delim_whitespace=True).iloc[:,1:]
        return df

    if residues is None:
        topology = load_topology(top)
        residues = [i+1 for i in range(1, topology.n_residues-1)]

    with open(input_file, 'w') as file:
        file.writelines([f'parm {top}\n',
                         f'trajin {traj}\n'])
        for res in residues:
            file.writelines([f'dihedral {res}_phi :{res-1}@C :{res}@N :{res}@CA :{res}@C out {output_file}\n',
                             f'dihedral {res}_psi :{res}@N :{res}@CA :{res}@C :{res+1}@N out {output_file}\n',
                             f'dihedral {res}_chi1 :{res}@N :{res}@CA :{res}@CB :{res}@CG out {output_file}\n']) # for residues that do not have @CG, this will just be ignored
        file.writelines(['go\nquit\n'])

    run(['cpptraj', '-i', input_file], stdout=DEVNULL)

    df = read_csv(output_file, delim_whitespace=True).iloc[:,1:]

    return df


def __filter(traj, ref, divergences, filter_val, filter_type):
    # return everything if no filter
    if filter_val is None:
        return traj, ref, divergences
    
    # if based on cutoff
    # get number of values that match
    if filter_type == 'cutoff':
        filter_val = len(where(divergences>=filter_val)[0])
    # if get top values and fraction
    # convert to number
    elif filter_type == 'top' and isinstance(filter_val, float):
        filter_val = int(len(divergences)*filter_val)
    # finally get indexes that are sorted and match criteria
    idxs = divergences.argsort()[::-1][:filter_val]
    # reference might not have all the same masks as trajectory
    ref_idxs = array([where(ref.columns == traj.columns[idx])[0][0] for idx in idxs if traj.columns[idx] in ref.columns])

    return traj.iloc[:,idxs], ref.iloc[:,ref_idxs], divergences[idxs]
    

def __save_frame(trajectory, topology, workdir):
    with open(f'{workdir}/save.in', 'w') as file:
        file.writelines([f'parm {topology}\n',
                         f'trajin {trajectory} 1 1\n',
                         f'trajout {workdir}/traj_f1.pdb\n',
                          'go\nquit\n'])
    
    run(['cpptraj', '-i', f'{workdir}/save.in'], stdout=DEVNULL)

    os.remove(f'{workdir}/save.in')

    return f'{workdir}/traj_f1.pdb'


def __plot_distributions(traj_distances, ref_distances, hbond_divergences, traj_dihedrals, ref_dihedrals, dihedral_divergences, colors, plot):
    n_cols = ((len(hbond_divergences) + len(dihedral_divergences))//3)+1
    height = n_cols*5
    fig, ax = _subplots(n_cols, 3, figsize=(15, height))
    ax = hstack(ax)
    for i in range(len(ax)):
        if i < len(hbond_divergences):
            _kdeplot(traj_distances.iloc[:,i], ax=ax[i], shade=True, color=colors[0], label='trajectory')
            if traj_distances.columns[i] in ref_distances.columns:
                idx = where(ref_distances.columns==traj_distances.columns[i])[0][0]
                _kdeplot(ref_distances.iloc[:,idx], ax=ax[i], shade=True, color=colors[1], label='reference')
            ax[i].set_xlabel('distance/$\AA$', size=14)
            ax[i].set_ylabel('density', size=14)
            ax[i].set_title(f'{traj_distances.columns[i]} (freq diff={round(hbond_divergences[i], 2)}%)')
            ax[i].legend()
        elif i-len(hbond_divergences) < len(dihedral_divergences):
            j = i-len(hbond_divergences)
            _kdeplot(traj_dihedrals.iloc[:,j], ax=ax[i], shade=True, color=colors[0], label='trajectory')
            _kdeplot(ref_dihedrals.iloc[:,j], ax=ax[i], shade=True, color=colors[1], label='reference')
            ax[i].set_xlabel('dihedral/$^{\circ}$', size=14)
            ax[i].set_ylabel('density', size=14)
            ax[i].set_title(f'{traj_dihedrals.columns[j]} (div={round(dihedral_divergences[j], 2)})')
            ax[i].legend()
        else:
            ax[i].set_visible(False)
    fig.tight_layout()

    fig.savefig(plot, dpi=300)
    
    return fig, ax


def __generate_pymol(hbonds, hbond_differences, dihedrals, pdb, workdir, fname):
    commands = [f'load {pdb}, frame1\n',
                 'color gray80, elem C\n',
                 'bg_color white\n']

    # show Hbonds
    for i in range(len(hbonds)):
        # display residues
        atom_masks = hbonds[i].split('-')
        residues = [mask.split('_')[1].split('@')[0] for mask in atom_masks[:2]]
        commands += [f'show sticks, resi {res}\ncolor yellow, resi {res} and elem C\n' for res in residues]
        # show bond
        atom1 = atom_masks[0].split('_')[1].split('@')[1]
        atom2 = atom_masks[2]
        commands += [f'distance d{i+1}, {residues[0]}/{atom1}, {residues[1]}/{atom2}\n']
        # color bond based on whether it is more frequent in trajectory or reference
        if hbond_differences[i] <= 0:
            color = 'red'
        else:
            color = 'green'
        commands += [f'hide label, d{i+1}\n', f'color {color}, d{i+1}\n']
    
    # show dihedrals
    for i in range(len(dihedrals)):
        # show residue
        res = int(dihedrals[i].split('_')[0])
        commands += [f'show sticks, resi {res}\ncolor purple, resi {res} and elem C\n']
        # show dihedral
        dih_type = dihedrals[i].split('_')[1]
        if dih_type == 'phi':
            mask = f'{res-1}/C, {res}/N, {res}/CA, {res}/C'
        elif dih_type == 'psi':
            mask = f'{res}/N, {res}/CA, {res}/C, {res+1}/N'
        else:
            mask = f'{res}/N, {res}/CA, {res}/CB, {res}/CG'
        commands += [f'dihedral t{i+1}, {mask}\n', f'hide label, t{i+1}\n']

    commands += [f'save {workdir}/{fname}.pse\n']

    with open(f'{workdir}/generate_session.pml', 'w') as file:
        file.writelines(commands)
    run(['pymol', '-c', f'{workdir}/generate_session.pml'], stdout=DEVNULL)

    return commands
        