import numpy as np
from mods import physics
from mods.ecal import main as emain
cellMap = np.loadtxt('./mods/ecal/cellmodule.txt') # Add bdt/ for running batch

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

'''
def recoil_tracks(f_dict, args, e_store, lq):

    """ Number of tracks and their eneries """
    # Consider minimu pmag to count

    n_recoil_tracks = 0
    max_pmag = 0
    for ecalSPHit in next( iter( args.values() ) ):

            if ecalSPHit.getCharge() == 0: continue

            n_recoil_tracks += 1
            if physics.mag( ecalSPHit.getMomentum() ) > max_pmag:
                max_pmag = physics.mag( ecalSPHit.getMomentum() )

    f_dict['n_recoil_tracks'] = n_recoil_tracks
    f_dict['max_track_pmag'] = max_pmag
'''

def recoilPT(f_dict, args, e_store, lq):

    """ Save recoilPMag from ecal_init to tree """

    f_dict['recoilPT'] = e_store['recoilPMag']

def decay_verts(f_dict, args, x_store, lq):

    """
    A' decay verticies (only applicable to vissig samples)
    +/- ~1 mm as really using electron creation verticies
    """

    Abi = 1
    for simParticle in next( iter( args.values() ) ):
        Abi += 1
        if Abi != 4: continue
        #f_dict['Ap_dt'] = simParticle[1].getTime Is overloading not sig anw
        f_dict['Ap_dx'], f_dict['Ap_dy'], f_dict['Ap_dz'] \
                = simParticle[1].getVertex()
