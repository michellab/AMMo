import pytraj as pt


def __validate_input(trajectory, topology):
    if isinstance(trajectory, str):
        if topology is not None:
            trajectory = pt.load(trajectory, top=topology)
        else:
            raise ValueError('if trajectory is "str", topology must be provided')
    elif not isinstance(trajectory, pt.Trajectory):
        raise TypeError('trajectory has to be of type "str" or "pytraj.Trajectory"')
        
    return trajectory


def __rmsd(trajectory, mask, reference, shared):
    if reference is None:
        raise ValueError('Please provide a reference for RMSD calculation')
    else:
        reference = pt.load(reference)
        
    if shared is None:
        shared = '*'
    
    aligned = pt.align(trajectory, shared, ref=reference, ref_mask=shared)
    
    rmsd = pt.rmsd(aligned, mask, ref=reference, ref_mask=mask, nofit=True, update_coordinate=False)
    return rmsd
        

def __distance(trajectory, mask):
    return pt.calc_distance(trajectory, mask)


def __torsion(trajectory, mask):
    return pt.calc_dihedral(trajectory, mask)


def featurize(trajectory, feature, mask, reference=None, shared='!@/H', topology=None):
    """Reduce trajectory data to a single feature.
    
    Parameters
    ----------
    trajectory : str, pytraj.Trajectory
        trajectory to be featurized. If 'str', 'topology' must be specified
        and it will be loaded. 
       
    feature : str
        feature type. Allowed: 'rmsd', 'distance', 'torsion'
        
    mask : str
        AMBER selection mask for the feature
        
    reference : str, pytraj.Trajectory
        reference for RMSD calculations. If 'str' of a PDB file, it will be
        loaded
        
    shared : str
        AMBER selection mask for shared atoms between trajectory and reference.
        If 'None', all atoms will be used. Default is not H, meaning all heavy atoms
        
    topology : str
        topology of the system
        
    Returns
    -------
    featurized : numpy.array
        a time series of the feature
    """
    trajectory = __validate_input(trajectory, topology)
        
    # calculate the feature
    if feature == 'rmsd':
        featurized = __rmsd(trajectory, mask, reference, shared)
    elif feature == 'distance':
        featurized = __distance(trajectory, mask)
    elif feature == 'torsion':
        featurized = __torsion(trajectory, mask)
    else:
        raise ValueError('"feature" must be "rmsd", "distance" or "torsion"')
        
    return featurized