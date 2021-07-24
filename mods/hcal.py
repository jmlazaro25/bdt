import math
import numpy as np

# Geometry Constants
# Layer 1 has no absorber, layers 2 and 3 have absorber of different thickness
air_thick = 2
scint_thick = 20
back_dx = 3100
back_dy = 3100
back_numLayers1 = 0
back_numLayers2 = 100
back_numLayers3 = 0
back_abso2_thick = 25
back_abso3_thick = 50
back_layer1_thick = scint_thick + air_thick
back_layer2_thick = back_abso2_thick + scint_thick + 2.0*air_thick
back_layer3_thick = back_abso3_thick + scint_thick + 2.0*air_thick
back_dz1 = back_numLayers1*back_layer1_thick
back_dz2 = back_numLayers2*back_layer2_thick
back_dz3 = back_numLayers3*back_layer3_thick
hcal_back_dz = back_dz1 + back_dz2 + back_dz3

# Side HCal Layer component
sideTB_layers = 28
sideLR_layers = 26
side_abso_thick = 20
hcal_side_dz = 600

# Macro
back_start_z = 860.5 # ecal_front_z + hcal_side_dz + 20
hcal_dz = hcal_back_dz + hcal_side_dz

# DetDescr
# HcalSection BACK = 0, TOP = 1, BOTTOM = 2, LEFT = 4, RIGHT = 3
SECTION_MASK = 0x7 # space for up to 7 sections
SECTION_SHIFT = 18
LAYER_MASK = 0xFF  # space for up to 255 layers
LAYER_SHIFT = 10
STRIP_MASK = 0xFF  # space for 255 strips/layer
STRIP_SHIFT = 0

# Analysis Constants
maxTime = 50
maxDepth = 4000
back_segments = 3

##################################################
# Global functions
##################################################

# ID-related info
##################################################

def section(hit):

    """ Get sectionID from hcal hit """

    return (hit.getID() >> SECTION_SHIFT) & SECTION_MASK

def layer(hit):

    """ Get layerID from hcal hit """

    return (hit.getID() >> LAYER_SHIFT) & LAYER_MASK

def strip(hit):

    """ Get stripID from hcal hit """

    return (hit.getID() >> STRIP_SHIFT) & STRIP_MASK

# Global features
##################################################

def maxPE(f_dict, args, h_store, lq):

    """ Copy maxPE from LDMX_events """

    hcalVeto = next( iter( args.values() ) )
    hcalVeto.getMaxPEHit().getPE()

    if 0 <= hcalVeto.getMaxPEHit().getPE() < 1e4:
        f_dict['maxPE'] = hcalVeto.getMaxPEHit().getPE()
        h_store['searchMaxPE'] = False

    else: h_store['searchMaxPE'] = True # maxPE left at default 0 at this stage

def maxPEsearch(f_dict, args, h_store, hcalRecHit):

    """ If needed, search for maxPE """

    if h_store['searchMaxPE']:

        if f_dict['maxPE'] < hcalRecHit.getPE() < 1e4:
            f_dict['maxPE'] = hcalRecHit.getPE()

# Consider instead section arg for all functions
##################################################
# Back functions
##################################################

def prep_hcal_lfs(f_dict, args, h_store, lq):

    """ Just add back_hits to store (othewise it'll keep being overwritten) """

    h_store['back_hits'] = np.zeros(5)

def collect(f_dict, args, h_store, hcalRecHit):

    """ Collect hcal hit info in desired array format """

    # Hit and energy-weighted hit (info) arrays
    if not section(hcalRecHit) == 0: return
    if not hcalRecHit.getEnergy() > 0: return
    if hcalRecHit.getTime() >= maxTime: return
    if hcalRecHit.getZPos() > maxDepth: return

    hcalRecHit = np.array( [
                    hcalRecHit.getEnergy(),
                    hcalRecHit.getXPos(),
                    hcalRecHit.getYPos(),
                    hcalRecHit.getZPos(),
                    hcalRecHit.getPE()
                    ] )

    h_store['back_hits'] = np.vstack( (h_store['back_hits'], hcalRecHit) )

def back_v1_all(f_dict, args, h_store, lq):

    """ Calculate collective part of 'first' set of hcal features  """
    
    # Don't try if there are no hits
    if type(h_store['back_hits'][0]) != np.ndarray:
        return

    f_dict['back_nHits'] = len(h_store['back_hits']) - 1
    f_dict['back_totE'] = sum( h_store['back_hits'][1:,0] )
    f_dict['back_totPE'] = sum( h_store['back_hits'][1:,4] )
    f_dict['back_maxE'] = max( h_store['back_hits'][1:,0] )
    f_dict['back_maxPE'] =  max( h_store['back_hits'][1:,4] )
    f_dict['back_avgE'] = np.mean( h_store['back_hits'][1:,0] )
    h_store['back_avg_e_x'] = np.average(
                            h_store['back_hits'][1:,1],
                            weights = h_store['back_hits'][1:,0]
                            )
    h_store['back_avg_e_y'] = np.average(
                            h_store['back_hits'][1:,2],
                            weights = h_store['back_hits'][1:,0]
                            )
    h_store['back_avg_e_z'] = np.average(
                            h_store['back_hits'][1:,3], 
                            weights = h_store['back_hits'][1:,0]
                            )
    f_dict['back_avgPE'] = np.mean( h_store['back_hits'][1:,4] )
    f_dict['back_avgR'] = np.average(
            np.sqrt(
                (h_store['back_hits'][1:,1]-h_store['back_avg_e_x'])**2 +\
                (h_store['back_hits'][1:,2]-h_store['back_avg_e_y'])**2
                ),
            weights = h_store['back_hits'][1:,0]
            )
    f_dict['back_std_e_x'] = math.sqrt(
            np.average( 
                (h_store['back_hits'][1:,1] - h_store['back_avg_e_x'])**2 ,
                weights = h_store['back_hits'][1:,0]
                )
            )
    f_dict['back_std_e_y'] = math.sqrt( 
            np.average(
                (h_store['back_hits'][1:,2] - h_store['back_avg_e_y'])**2 ,
                weights = h_store['back_hits'][1:,0]
                )
            )
    f_dict['back_std_e_z'] = math.sqrt(
            np.average(
                (h_store['back_hits'][1:,3] - h_store['back_avg_e_z'])**2 ,
                weights = h_store['back_hits'][1:,0] 
                )
            )
    f_dict['back_dz_e'] = max( h_store['back_hits'][1:,3] )\
                            - h_store['back_avg_e_z']

def back_v1_seg(f_dict, args, h_store, lq):

    """ Calculate regional parts of 'first' set of hcal features """

    # Don't try if there are no hits
    if type(h_store['back_hits'][0]) != np.ndarray:
        return

     # Regional Stats
    hits_e_1 = hits_e_2 = hits_e_3 = np.zeros(5)
    
    # Group hits into regions
    for hcalRecHit in h_store['back_hits']:
        if hcalRecHit[3] < h_store['back_avg_e_z']:
            hits_e_1 = np.vstack( (hits_e_1, hcalRecHit) )
        elif hcalRecHit[3] < h_store['back_avg_e_z'] + f_dict['back_std_e_z']:
            hits_e_2 = np.vstack( (hits_e_2, hcalRecHit) )
        else:
            hits_e_3 = np.vstack( (hits_e_3, hcalRecHit) )

    # For each region, find averages, radii, etc.
    regions = [0, hits_e_1, hits_e_2, hits_e_3]
    for i in range(1,4):
        if type(regions[i][0]) == np.ndarray: # Make sure the region has any hits
            f_dict['back_nHits_{}e'.format(i)] = len( regions[i] ) - 1
            f_dict['back_tot_{}e_e'.format(i)] = sum( regions[i][1:,0] )
            if f_dict['back_tot_{}e_e'.format(i)] == 0: continue # Rounding round error
            f_dict['back_tot_{}e_pe'.format(i)] = sum( regions[i][1:,4] )
            f_dict['back_max_{}e_e'.format(i)] = max( regions[i][1:,0] )
            f_dict['back_max_{}e_pe'.format(i)] = max( regions[i][1:,4] )
            f_dict['back_avg_{}e_e'.format(i)] = np.mean( regions[i][1:,0] )
            avg_e_x = np.average( regions[i][1:,1], weights = regions[i][1:,0] )
            avg_e_y = np.average( regions[i][1:,2], weights = regions[i][1:,0] )
            f_dict['back_avg_{}e_pe'.format(i)] = np.mean( regions[i][1:,0] )
            f_dict['back_avg_{}e_r'.format(i)] = np.average(
                    np.sqrt(
                        (regions[i][1:,1]-avg_e_x)**2 +\
                        (regions[i][1:,2]-avg_e_y)**2
                        ),
                    weights = regions[i][1:,0]
                    )
            if f_dict['back_nHits_{}e'.format(i)] > 1: # Std = 0 if there's only 1 hit
                f_dict['back_std_{}e_e'.format(i)] = math.sqrt(
                        sum((regions[i][1:,0] - f_dict['back_avg_{}e_e'.format(i)])**2)/\
                            f_dict['back_nHits_{}e'.format(i)] )
                f_dict['back_std_{}e_x'.format(i)] = math.sqrt( np.average( 
                    (regions[i][1:,1] - avg_e_x)**2,
                    weights = regions[i][1:,0] ) )
                f_dict['back_std_{}e_y'.format(i)] = math.sqrt( np.average( 
                    (regions[i][1:,2] - avg_e_y)**2,
                    weights = regions[i][1:,0] ) )
                f_dict['back_std_{}e_pe'.format(i)] = math.sqrt(
                        sum((regions[i][1:,4] - f_dict['back_avg_{}e_pe'.format(i)])**2)/\
                            f_dict['back_nHits_{}e'.format(i)] )

##################################################
# Side functions
##################################################
