import math
import numpy as np
from mods import physics
import mods.ecal.main as emain

# Constants
segLayers = [0, 6, 17, 34]
nSegments = len(segLayers) - 1

# Segment-containment feature functions

def segcont_init(f_dict, args, e_store, lq):

    """ Use recoilPMag and recoilTheta to select containment radii """

    # Set electron RoC binnings
    if e_store['globals']['ebeam'] == 4:

        e_store['e_radii'] = emain.radius68_thetalt10_plt500
        if e_store['recoilTheta'] < 10 and e_store['recoilPMag'] >= 500:
            e_store['e_radii'] = emain.radius68_thetalt10_pgt500
        elif 10 <= e_store['recoilTheta'] < 20:
            e_store['e_radii'] = emain.radius68_theta10to20
        elif e_store['recoilTheta'] >= 20:
            e_store['e_radii'] = emain.radius68_thetagt20

        # Always use default binning for photon RoC
        e_store['g_radii'] = emain.radius68_thetalt10_plt500
        return

    # If here, assuming 8 GeV
    e_store['e_radii'] = emain.radius68_thetalt6_plt1000
    if e_store['recoilTheta'] < 6 and e_store['recoilPMag'] >= 1000:
        e_store['e_radii'] = emain.radius68_thetalt6_pgt1000
    elif 6 <= e_store['recoilTheta'] < 20:
        e_store['e_radii'] = emain.radius68_theta6to15
    elif e_store['recoilTheta'] >= 20:
        e_store['e_radii'] = emain.radius68_thetagt20

    e_store['g_radii'] = emain.radius68_thetalt6_plt1000

def prep_ecal_lfs(f_dict, args, e_store, lq):

    """ Just add back_hits to store (othewise it'll keep being overwritten) """

    e_store['ecal_hits'] = np.zeros(4)

def collect(f_dict, args, e_store, ecalRecHit):

    """ Collect ecal hit info in desired array format """

    # Hit and energy-weighted hit (info) arrays
    if not ecalRecHit.getEnergy() > 0: return

    ecalRecHit = np.array( [
                    ecalRecHit.getEnergy(),
                    ecalRecHit.getXPos(),
                    ecalRecHit.getYPos(),
                    emain.layer(ecalRecHit)
                    ] )

    e_store['ecal_hits'] = np.vstack( (e_store['ecal_hits'], ecalRecHit) )

def rsegcont(f_dict, args, e_store, lq):

    """ Calculate relative segment features """

    # Don't try if there are no hits
    if type(e_store['ecal_hits'][0]) != np.ndarray: return

    # Regional stats go here
    hits1 = hits2 = hits3 = np.zeros(4)

    # Group hits into regions
    for ecalRecHit in e_store['ecal_hits'][1:]:
        if ecalRecHit[3] < e_store['avgLayerHit']:
            hits1 = np.vstack( (hits1, ecalRecHit) )
        elif ecalRecHit[3] < e_store['avgLayerHit'] + f_dict['stdLayerHit']:
            hits2 = np.vstack( (hits2, ecalRecHit) )
        else:
            hits3 = np.vstack( (hits3, ecalRecHit) )

    # For each region, find averages, radii, etc.
    regions = [0, hits1, hits2, hits3]
    for s in range(1,4):
        if type(regions[s][0]) == np.ndarray: # Make sure the region has any hits
            f_dict[f'nHits_rs{s}'] = len( regions[s] ) - 1
            f_dict[f'energy_rs{s}'] = sum( regions[s][1:,0] )
            #if f_dict['back_tot_{}e_e'.format(i)] == 0: continue # Rounding round error
            """
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
            """

def segcont_means(f_dict, args, e_store, ecalRecHit):

    """ Segment-containment calculations for first loop over EcalRecHits """

    if ecalRecHit.getEnergy() <= 0: return

    # Reused often
    energy = ecalRecHit.getEnergy()
    layer = emain.layer(ecalRecHit)
    xy_pair = ( ecalRecHit.getXPos(), ecalRecHit.getYPos() )

    # Distance to electron trajectory
    if e_store['e_traj'] != None:
        xy_e_traj = (e_store['e_traj'][layer][0], e_store['e_traj'][layer][1])
        distance_e_traj = physics.dist(xy_pair, xy_e_traj)
    else: distance_e_traj = -1.0

    # Distance to photon trajectory
    if e_store['g_traj'] != None:
        xy_g_traj = (e_store['g_traj'][layer][0], e_store['g_traj'][layer][1])
        distance_g_traj = physics.dist(xy_pair, xy_g_traj)
    else: distance_g_traj = -1.0

    # Decide which longitudinal segment the hit is in and add to sums
    for s in range(1, nSegments + 1):

        if not ( segLayers[s - 1] <= layer and (layer <= segLayers[s] - 1) ):
            continue

        f_dict['energy_s{}'.format(s)] += energy
        f_dict['nHits_s{}'.format(s)] += 1
        f_dict['xMean_s{}'.format(s)] += xy_pair[0]*energy
        f_dict['yMean_s{}'.format(s)] += xy_pair[1]*energy
        f_dict['layerMean_s{}'.format(s)] += layer*energy

        # Decide which containment region the hit is in and add to sums
        for r in range(1, emain.contRegions + 1):

            if (r - 1)*e_store['e_radii'][layer] <= distance_e_traj \
                    and distance_e_traj < r*e_store['e_radii'][layer]:
                f_dict['eContEnergy_x{}_s{}'.format(r,s)] += energy
                f_dict['eContNHits_x{}_s{}'.format(r,s)] += 1
                f_dict['eContXMean_x{}_s{}'.format(r,s)] += xy_pair[0]*energy
                f_dict['eContYMean_x{}_s{}'.format(r,s)] += xy_pair[1]*energy
                f_dict['eContLayerMean_x{}_s{}'.format(r,s)] += layer*energy

            if (r - 1)*e_store['g_radii'][layer] <= distance_g_traj \
                    and distance_g_traj < r*e_store['g_radii'][layer]:
                f_dict['gContEnergy_x{}_s{}'.format(r,s)] += energy
                f_dict['gContNHits_x{}_s{}'.format(r,s)] += 1
                f_dict['gContXMean_x{}_s{}'.format(r,s)] += xy_pair[0]*energy
                f_dict['gContYMean_x{}_s{}'.format(r,s)] += xy_pair[1]*energy
                f_dict['gContLayerMean_x{}_s{}'.format(r,s)] += layer*energy

            if distance_e_traj > r*e_store['e_radii'][layer] \
                    and distance_g_traj > r*e_store['g_radii'][layer]:
                f_dict['oContEnergy_x{}_s{}'.format(r,s)] += energy
                f_dict['oContNHits_x{}_s{}'.format(r,s)] += 1
                f_dict['oContXMean_x{}_s{}'.format(r,s)] += xy_pair[0]*energy
                f_dict['oContYMean_x{}_s{}'.format(r,s)] += xy_pair[1]*energy
                f_dict['oContLayerMean_x{}_s{}'.format(r,s)] += layer*energy

def segcont_means_norm(f_dict, args, e_store, lq):

    """ Normalize mean calculations started in segcont_means """

    for s in range(1, nSegments + 1):

        if f_dict['energy_s{}'.format(s)] > 0:
            f_dict['xMean_s{}'.format(s)] /= f_dict['energy_s{}'.format(s)]
            f_dict['yMean_s{}'.format(s)] /= f_dict['energy_s{}'.format(s)]
            f_dict['layerMean_s{}'.format(s)] /= f_dict['energy_s{}'.format(s)]

        for r in range(1, emain.contRegions + 1):

            if f_dict['eContEnergy_x{}_s{}'.format(r,s)] > 0:
                f_dict['eContXMean_x{}_s{}'.format(r,s)] \
                        /= f_dict['eContEnergy_x{}_s{}'.format(r,s)]
                f_dict['eContYMean_x{}_s{}'.format(r,s)] \
                        /= f_dict['eContEnergy_x{}_s{}'.format(r,s)]
                f_dict['eContLayerMean_x{}_s{}'.format(r,s)] \
                        /= f_dict['eContEnergy_x{}_s{}'.format(r,s)]

            if f_dict['gContEnergy_x{}_s{}'.format(r,s)] > 0:
                f_dict['gContXMean_x{}_s{}'.format(r,s)] \
                        /= f_dict['gContEnergy_x{}_s{}'.format(r,s)]
                f_dict['gContYMean_x{}_s{}'.format(r,s)] \
                        /= f_dict['gContEnergy_x{}_s{}'.format(r,s)]
                f_dict['gContLayerMean_x{}_s{}'.format(r,s)] \
                        /= f_dict['gContEnergy_x{}_s{}'.format(r,s)]

            if f_dict['oContEnergy_x{}_s{}'.format(r,s)] > 0:
                f_dict['oContXMean_x{}_s{}'.format(r,s)] \
                        /= f_dict['oContEnergy_x{}_s{}'.format(r,s)]
                f_dict['oContYMean_x{}_s{}'.format(r,s)] \
                        /= f_dict['oContEnergy_x{}_s{}'.format(r,s)]
                f_dict['oContLayerMean_x{}_s{}'.format(r,s)] \
                        /= f_dict['oContEnergy_x{}_s{}'.format(r,s)]

def segcont_stds(f_dict, args, e_store, ecalRecHit):

    """ Segment-containment calculations for second loop over EcalRecHits """

    if ecalRecHit.getEnergy() <= 0: return

    # Reused often
    energy = ecalRecHit.getEnergy()
    layer = emain.layer(ecalRecHit)
    xy_pair = ( ecalRecHit.getXPos(), ecalRecHit.getYPos() )

    # Distance to electron trajectory
    if e_store['e_traj'] != None:
        xy_e_traj = (e_store['e_traj'][layer][0], e_store['e_traj'][layer][1])
        distance_e_traj = physics.dist(xy_pair, xy_e_traj)
    else: distance_e_traj = -1.0

    # Distance to photon trajectory
    if e_store['g_traj'] != None:
        xy_g_traj = (e_store['g_traj'][layer][0], e_store['g_traj'][layer][1])
        distance_g_traj = physics.dist(xy_pair, xy_g_traj)
    else: distance_g_traj = -1.0

    # Decide which longitudinal segment the hit is in and add to sums
    for s in range(1, nSegments + 1):

        if not ( segLayers[s - 1] <= layer and (layer <= segLayers[s] - 1) ):
            continue

        f_dict['xStd_s{}'.format(s)] += ( (xy_pair[0] \
                - f_dict['xMean_s{}'.format(s)])**2 )*energy
        f_dict['yStd_s{}'.format(s)] += ( (xy_pair[1] \
                - f_dict['yMean_s{}'.format(s)])**2 )*energy
        f_dict['layerStd_s{}'.format(s)] += ( (layer \
                - f_dict['layerMean_s{}'.format(s)])**2 )*energy

        # Decide which containment region the hit is in and add to sums
        for r in range(1, emain.contRegions + 1):

            if (r - 1)*e_store['e_radii'][layer] <= distance_e_traj \
              and distance_e_traj < r*e_store['e_radii'][layer]:
                f_dict['eContXStd_x{}_s{}'.format(r,s)] += ((xy_pair[0] \
                        - f_dict['eContXMean_x{}_s{}'.format(r,s)])**2)*energy
                f_dict['eContYStd_x{}_s{}'.format(r,s)] += ((xy_pair[1] \
                        - f_dict['eContYMean_x{}_s{}'.format(r,s)])**2)*energy
                f_dict['eContLayerStd_x{}_s{}'.format(r,s)] += ((layer \
                        - f_dict['eContLayerMean_x{}_s{}'.format(r,s)])**2)\
                        *energy

            if (r - 1)*e_store['g_radii'][layer] <= distance_g_traj \
              and distance_g_traj < r*e_store['g_radii'][layer]:
                f_dict['gContXStd_x{}_s{}'.format(r,s)] += ((xy_pair[0] \
                        - f_dict['gContXMean_x{}_s{}'.format(r,s)])**2)*energy
                f_dict['gContYStd_x{}_s{}'.format(r,s)] += ((xy_pair[1] \
                        - f_dict['gContYMean_x{}_s{}'.format(r,s)])**2)*energy
                f_dict['gContLayerStd_x{}_s{}'.format(r,s)] += ((layer \
                        - f_dict['gContLayerMean_x{}_s{}'.format(r,s)])**2)\
                        *energy

            if distance_e_traj > r*e_store['e_radii'][layer] \
              and distance_g_traj > r*e_store['g_radii'][layer]:
                f_dict['oContXStd_x{}_s{}'.format(r,s)] += ((xy_pair[0] \
                        - f_dict['oContXMean_x{}_s{}'.format(r,s)])**2)*energy
                f_dict['oContYStd_x{}_s{}'.format(r,s)] += ((xy_pair[1] \
                        - f_dict['oContYMean_x{}_s{}'.format(r,s)])**2)*energy
                f_dict['oContLayerStd_x{}_s{}'.format(r,s)] += ((layer \
                        - f_dict['oContLayerMean_x{}_s{}'.format(r,s)])**2)\
                        *energy

def segcont_stds_norm(f_dict, args, e_store, lq):

    """ Normalize standard deviation calculations started in seccont_stds """

    for s in range(1, nSegments + 1):

        if f_dict['energy_s{}'.format(s)] > 0:
            f_dict['xStd_s{}'.format(s)] = math.sqrt(
                    f_dict['xStd_s{}'.format(s)] \
                    / f_dict['energy_s{}'.format(s)]    )
            f_dict['yStd_s{}'.format(s)] = math.sqrt(
                    f_dict['yStd_s{}'.format(s)] \
                    / f_dict['energy_s{}'.format(s)]    )
            f_dict['layerStd_s{}'.format(s)] = math.sqrt(
                    f_dict['layerStd_s{}'.format(s)] \
                    / f_dict['energy_s{}'.format(s)]    )

        for r in range(1, emain.contRegions + 1):

            if f_dict['eContEnergy_x{}_s{}'.format(r,s)] > 0:
                f_dict['eContXStd_x{}_s{}'.format(r,s)] \
                        = math.sqrt(f_dict['eContXStd_x{}_s{}'.format(r,s)] \
                        / f_dict['eContEnergy_x{}_s{}'.format(r,s)] )
                f_dict['eContYStd_x{}_s{}'.format(r,s)] \
                        = math.sqrt(f_dict['eContYStd_x{}_s{}'.format(r,s)] \
                        / f_dict['eContEnergy_x{}_s{}'.format(r,s)] )
                f_dict['eContLayerStd_x{}_s{}'.format(r,s)] \
                        = math.sqrt(f_dict['eContLayerStd_x{}_s{}'\
                                                                .format(r,s)] \
                                / f_dict['eContEnergy_x{}_s{}'.format(r,s)] )

            if f_dict['gContEnergy_x{}_s{}'.format(r,s)] > 0:
                f_dict['gContXStd_x{}_s{}'.format(r,s)] \
                        = math.sqrt(f_dict['gContXStd_x{}_s{}'.format(r,s)] \
                        / f_dict['gContEnergy_x{}_s{}'.format(r,s)] )
                f_dict['gContYStd_x{}_s{}'.format(r,s)] \
                        = math.sqrt(f_dict['gContYStd_x{}_s{}'.format(r,s)] \
                        / f_dict['gContEnergy_x{}_s{}'.format(r,s)] )
                f_dict['gContLayerStd_x{}_s{}'.format(r,s)] \
                        = math.sqrt(f_dict['gContLayerStd_x{}_s{}'\
                                                                .format(r,s)] \
                                / f_dict['gContEnergy_x{}_s{}'.format(r,s)] )

            if f_dict['oContEnergy_x{}_s{}'.format(r,s)] > 0:
                f_dict['oContXStd_x{}_s{}'.format(r,s)] \
                        = math.sqrt(f_dict['oContXStd_x{}_s{}'.format(r,s)] \
                        / f_dict['oContEnergy_x{}_s{}'.format(r,s)] )
                f_dict['oContYStd_x{}_s{}'.format(r,s)] \
                        = math.sqrt(f_dict['oContYStd_x{}_s{}'.format(r,s)] \
                        / f_dict['oContEnergy_x{}_s{}'.format(r,s)] )
                f_dict['oContLayerStd_x{}_s{}'.format(r,s)] \
                        = math.sqrt(f_dict['oContLayerStd_x{}_s{}'\
                                                                .format(r,s)] \
                                / f_dict['oContEnergy_x{}_s{}'.format(r,s)] )
