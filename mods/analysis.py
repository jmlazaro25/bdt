import numpy as np
from mods import physics
from mods.ecal import main as emain
cellMap = np.loadtxt('./mods/ecal/cellmodule.txt')

def fidcats(f_dict, args, e_store, lq):

    """ Determine if electron and/or photon is/are fiducial (bool) """

    if e_store['e_traj'] != None:
        for cell in cellMap:
            if physics.dist( cell[1:], e_store['e_traj'][0] ) <= emain.cell_radius:
                f_dict['e_fid'] = True
                break

    if e_store['g_traj'] != None:
        for cell in cellMap:
            if physics.dist( cell[1:], e_store['g_traj'][0] ) <= emain.cell_radius:
                f_dict['g_fid'] = True
                break

def trigger(f_dict, args, x_store, lq):

    """ Record if event cut by standard trigger """

    f_dict['trigPass'] = next( iter( args.values ) ).passed()

def recoilPT(f_dict, args, e_store, lq):

    """ Save recoilPMag from ecal_init to tree """

    f_dict['recoilPT'] = e_store['recoilPMag']
