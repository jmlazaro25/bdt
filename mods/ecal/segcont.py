import math
import numpy as np
import awkward as ak
from mods import physics
from copy import deepcopy
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

def rsegcont_means(f_dict, args, e_store, ecalRecHit):

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

        if s == 1 and not( layer < e_store['avgLayerHit'] ): continue
        elif s == 2 and not(
                            e_store['avgLayerHit'] \
                            < layer \
                            < e_store['avgLayerHit'] + f_dict['stdLayerHit']
                            ): continue
        elif s == 3 and not(
                            e_store['avgLayerHit'] + f_dict['stdLayerHit']
                            < layer): continue

        f_dict[f'energy_rs{s}'] += energy
        f_dict[f'nHits_rs{s}'] += 1
        f_dict[f'xMean_rs{s}'] += xy_pair[0]*energy
        f_dict[f'yMean_rs{s}'] += xy_pair[1]*energy
        e_store[f'layerMean_rs{s}'] += layer*energy

        # Decide which containment region the hit is in and add to sums
        for r in range(1, emain.contRegions + 1):

            if (r - 1)*e_store['e_radii'][layer] <= distance_e_traj \
                    and distance_e_traj < r*e_store['e_radii'][layer]:
                f_dict[f'eContEnergy_x{r}_rs{s}'] += energy
                f_dict[f'eContNHits_x{r}_rs{s}'] += 1
                f_dict[f'eContXMean_x{r}_rs{s}'] += xy_pair[0]*energy
                f_dict[f'eContYMean_x{r}_rs{s}'] += xy_pair[1]*energy
                e_store[f'eContLayerMean_x{r}_rs{s}'] += layer*energy

            if (r - 1)*e_store['g_radii'][layer] <= distance_g_traj \
                    and distance_g_traj < r*e_store['g_radii'][layer]:
                f_dict[f'gContEnergy_x{r}_rs{s}'] += energy
                f_dict[f'gContNHits_x{r}_rs{s}'] += 1
                f_dict[f'gContXMean_x{r}_rs{s}'] += xy_pair[0]*energy
                f_dict[f'gContYMean_x{r}_rs{s}'] += xy_pair[1]*energy
                e_store[f'gContLayerMean_x{r}_rs{s}'] += layer*energy

            if distance_e_traj > r*e_store['e_radii'][layer] \
                    and distance_g_traj > r*e_store['g_radii'][layer]:
                f_dict[f'oContEnergy_x{r}_rs{s}'] += energy
                f_dict[f'oContNHits_x{r}_rs{s}'] += 1
                f_dict[f'oContXMean_x{r}_rs{s}'] += xy_pair[0]*energy
                f_dict[f'oContYMean_x{r}_rs{s}'] += xy_pair[1]*energy
                e_store[f'oContLayerMean_x{r}_rs{s}'] += layer*energy

def rsegcont_means_norm(f_dict, args, e_store, lq):

    """ Normalize mean calculations started in segcont_means """

    for s in range(1, nSegments + 1):

        if f_dict[f'energy_rs{s}'] > 0:
            f_dict[f'xMean_rs{s}'] /= f_dict[f'energy_rs{s}']
            f_dict[f'yMean_rs{s}'] /= f_dict[f'energy_rs{s}']
            e_store[f'layerMean_rs{s}'] /= f_dict[f'energy_rs{s}']

        for r in range(1, emain.contRegions + 1):

            if f_dict[f'eContEnergy_x{r}_rs{s}'] > 0:
                f_dict[f'eContXMean_x{r}_rs{s}'] \
                        /= f_dict[f'eContEnergy_x{r}_rs{s}']
                f_dict[f'eContYMean_x{r}_rs{s}'] \
                        /= f_dict[f'eContEnergy_x{r}_rs{s}']
                e_store[f'eContLayerMean_x{r}_rs{s}'] \
                        /= f_dict[f'eContEnergy_x{r}_rs{s}']

            if f_dict[f'gContEnergy_x{r}_rs{s}'] > 0:
                f_dict[f'gContXMean_x{r}_rs{s}'] \
                        /= f_dict[f'gContEnergy_x{r}_rs{s}']
                f_dict[f'gContYMean_x{r}_rs{s}'] \
                        /= f_dict[f'gContEnergy_x{r}_rs{s}']
                e_store[f'gContLayerMean_x{r}_rs{s}'] \
                        /= f_dict[f'gContEnergy_x{r}_rs{s}']

            if f_dict[f'oContEnergy_x{r}_rs{s}'] > 0:
                f_dict[f'oContXMean_x{r}_rs{s}'] \
                        /= f_dict[f'oContEnergy_x{r}_rs{s}']
                f_dict[f'oContYMean_x{r}_rs{s}'] \
                        /= f_dict[f'oContEnergy_x{r}_rs{s}']
                e_store[f'oContLayerMean_x{r}_rs{s}'] \
                        /= f_dict[f'oContEnergy_x{r}_rs{s}']

def rsegcont_stds(f_dict, args, e_store, ecalRecHit):

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

        if s == 1 and not( layer < e_store['avgLayerHit'] ): continue
        elif s == 2 and not(
                            e_store['avgLayerHit'] \
                            < layer \
                            < e_store['avgLayerHit'] + f_dict['stdLayerHit']
                            ): continue
        elif s == 3 and not(
                            e_store['avgLayerHit'] + f_dict['stdLayerHit']
                            < layer): continue

        f_dict[f'xStd_rs{s}'] += ( (xy_pair[0] \
                - f_dict[f'xMean_rs{s}'])**2 )*energy
        f_dict[f'yStd_rs{s}'] += ( (xy_pair[1] \
                - f_dict[f'yMean_rs{s}'])**2 )*energy
        f_dict[f'layerStd_rs{s}'] += ( (layer \
                - e_store[f'layerMean_rs{s}'])**2 )*energy

        # Decide which containment region the hit is in and add to sums
        for r in range(1, emain.contRegions + 1):

            if (r - 1)*e_store['e_radii'][layer] <= distance_e_traj \
              and distance_e_traj < r*e_store['e_radii'][layer]:
                f_dict[f'eContXStd_x{r}_rs{s}'] += ((xy_pair[0] \
                        - f_dict[f'eContXMean_x{r}_rs{s}'])**2)*energy
                f_dict[f'eContYStd_x{r}_rs{s}'] += ((xy_pair[1] \
                        - f_dict[f'eContYMean_x{r}_rs{s}'])**2)*energy
                f_dict[f'eContLayerStd_x{r}_rs{s}'] += ((layer \
                        - e_store[f'eContLayerMean_x{r}_rs{s}'])**2)\
                        *energy

            if (r - 1)*e_store['g_radii'][layer] <= distance_g_traj \
              and distance_g_traj < r*e_store['g_radii'][layer]:
                f_dict[f'gContXStd_x{r}_rs{s}'] += ((xy_pair[0] \
                        - f_dict[f'gContXMean_x{r}_rs{s}'])**2)*energy
                f_dict[f'gContYStd_x{r}_rs{s}'] += ((xy_pair[1] \
                        - f_dict[f'gContYMean_x{r}_rs{s}'])**2)*energy
                f_dict[f'gContLayerStd_x{r}_rs{s}'] += ((layer \
                        - e_store[f'gContLayerMean_x{r}_rs{s}'])**2)\
                        *energy

            if distance_e_traj > r*e_store['e_radii'][layer] \
              and distance_g_traj > r*e_store['g_radii'][layer]:
                f_dict[f'oContXStd_x{r}_rs{s}'] += ((xy_pair[0] \
                        - f_dict[f'oContXMean_x{r}_rs{s}'])**2)*energy
                f_dict[f'oContYStd_x{r}_rs{s}'] += ((xy_pair[1] \
                        - f_dict[f'oContYMean_x{r}_rs{s}'])**2)*energy
                f_dict[f'oContLayerStd_x{r}_rs{s}'] += ((layer \
                        - e_store[f'oContLayerMean_x{r}_rs{s}'])**2)\
                        *energy

def rsegcont_stds_norm(f_dict, args, e_store, lq):

    """ Normalize standard deviation calculations started in seccont_stds """

    for s in range(1, nSegments + 1):

        if f_dict[f'energy_rs{s}'] > 0:
            f_dict[f'xStd_rs{s}'] = math.sqrt(
                    f_dict[f'xStd_rs{s}'] \
                    / f_dict[f'energy_rs{s}']    )
            f_dict[f'yStd_rs{s}'] = math.sqrt(
                    f_dict[f'yStd_rs{s}'] \
                    / f_dict[f'energy_rs{s}']    )
            f_dict[f'layerStd_rs{s}'] = math.sqrt(
                    f_dict[f'layerStd_rs{s}'] \
                    / f_dict[f'energy_rs{s}']    )

        for r in range(1, emain.contRegions + 1):

            if f_dict[f'eContEnergy_x{r}_rs{s}'] > 0:
                f_dict[f'eContXStd_x{r}_rs{s}'] \
                        = math.sqrt(f_dict[f'eContXStd_x{r}_rs{s}'] \
                        / f_dict[f'eContEnergy_x{r}_rs{s}'] )
                f_dict[f'eContYStd_x{r}_rs{s}'] \
                        = math.sqrt(f_dict[f'eContYStd_x{r}_rs{s}'] \
                        / f_dict[f'eContEnergy_x{r}_rs{s}'] )
                f_dict[f'eContLayerStd_x{r}_rs{s}'] \
                        = math.sqrt(f_dict[f'eContLayerStd_x{r}_rs{s}']\
                                / f_dict[f'eContEnergy_x{r}_rs{s}'] )

            if f_dict[f'gContEnergy_x{r}_rs{s}'] > 0:
                f_dict[f'gContXStd_x{r}_rs{s}'] \
                        = math.sqrt(f_dict[f'gContXStd_x{r}_rs{s}'] \
                        / f_dict[f'gContEnergy_x{r}_rs{s}'] )
                f_dict[f'gContYStd_x{r}_rs{s}'] \
                        = math.sqrt(f_dict[f'gContYStd_x{r}_rs{s}'] \
                        / f_dict[f'gContEnergy_x{r}_rs{s}'] )
                f_dict[f'gContLayerStd_x{r}_rs{s}'] \
                        = math.sqrt(f_dict[f'gContLayerStd_x{r}_rs{s}']\
                                / f_dict[f'gContEnergy_x{r}_rs{s}'] )

            if f_dict[f'oContEnergy_x{r}_rs{s}'] > 0:
                f_dict[f'oContXStd_x{r}_rs{s}'] \
                        = math.sqrt(f_dict[f'oContXStd_x{r}_rs{s}'] \
                        / f_dict[f'oContEnergy_x{r}_rs{s}'] )
                f_dict[f'oContYStd_x{r}_rs{s}'] \
                        = math.sqrt(f_dict[f'oContYStd_x{r}_rs{s}'] \
                        / f_dict[f'oContEnergy_x{r}_rs{s}'] )
                f_dict[f'oContLayerStd_x{r}_rs{s}'] \
                        = math.sqrt(f_dict[f'oContLayerStd_x{r}_rs{s}']\
                               / f_dict[f'oContEnergy_x{r}_rs{s}'] )

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

''' Massive irony...
or maybe not... (hard to tell speed of this vs method below without more tests
for which there is no time. This way at least it can't be any slower than 2X
dor doing segcont and rsegcont... right? (memory maketh :') )
def prep_ecal_lfs(f_dict, args, e_store, lq):

    """ Just add back_hits to store (othewise it'll keep being overwritten) """

    e_store['ecal_hits'] = []

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

    e_store['ecal_hits'].append( ecalRecHit )

def rsegcont(f_dict, args, e_store, lq):

    """ Calculate relative segment features """

    # Exclusive catagories in each containment region set
    cont = [[],[],[],[]] # [both, electron, gamma, outside]

    # Each rseg is a collection of multiples of cont (and there's 3 of them)
    rseg = [ deepcopy(cont) for _ in range(emain.contRegions) ]
    rsegs = [ deepcopy(rseg) for _ in range(nSegments) ]

    # Final struct note:
    # Relative_segments (3 for now)                                   (layer 1)
    # Each contain emain.ncontRegions (i.e. _x1, _x2, ..., _xn)       (layer 2)
    # Each _xi contains catagories (both, electron, gamma, outside)   (layer 3)
    # Each catagory contains hits (n of them)                         (layer 4)
    # Each hit is [energy, x, y, layer] (Ideally try z at some point) (layer 5)

    # Fill relative segments (and sub-arrays)
    for ecalRecHit in e_store['ecal_hits']:

        # Choose rsegment
        layer = int( ecalRecHit[3] )
        if layer < e_store['avgLayerHit']:
            rseg_i = 0
        elif layer < e_store['avgLayerHit'] + f_dict['stdLayerHit']:
            rseg_i = 1
        else:
            rseg_i = 2

        # Distance to electron trajectory
        if e_store['e_traj'] != None:
            xy_e_traj = (
                            e_store['e_traj'][layer][0],
                            e_store['e_traj'][layer][1]
                            )
            distance_e_traj = physics.dist( ecalRecHit[1:3] , xy_e_traj )
        else: distance_e_traj = -1.0

        # Distance to photon trajectory
        if e_store['g_traj'] != None:
            xy_g_traj = (
                            e_store['g_traj'][layer][0],
                            e_store['g_traj'][layer][1]
                            )
            distance_g_traj = physics.dist( ecalRecHit[1:3] , xy_g_traj )
        else: distance_g_traj = -1.0

        # Choose region(s)
        for r in range( emain.contRegions ):

            eRegion = gRegion = False

            if r*e_store['e_radii'][layer] <= distance_e_traj \
                    < (r + 1)*e_store['e_radii'][layer]:
                eRegion = True

            if r*e_store['g_radii'][layer] <= distance_g_traj \
                    < (r + 1)*e_store['g_radii'][layer]:
                gRegion = True

            # Do the actual filling
            if eRegion and gRegion: rsegs[rseg_i][r][0].append( ecalRecHit )
            elif eRegion: rsegs[rseg_i][r][1].append( ecalRecHit )
            elif gRegion: rsegs[rseg_i][r][2].append( ecalRecHit )
            else: rsegs[rseg_i][r][3].append( ecalRecHit )

    # For each section and region, find averages, stds, etc.
    # Probably should make a std function fo avoid the blocks below (lol)
    rsegs = ak.Array( rsegs )
    for s, rseg in enumerate(rsegs):

        # Overall segments Should probably remove x/yMeans for pT bias
        f_dict[f'nHits_rs{s + 1}'] = ak.count( rseg[0] ) // 4
        if f_dict[f'nHits_rs{s + 1}'] == 0: continue
        f_dict[f'energy_rs{s + 1}'] = ak.sum( rseg[0,:,:,0] )
        f_dict[f'xMean_rs{s + 1}'] = ak.mean( rseg[0,:,:,1], rseg[0,:,:,0] )
        f_dict[f'yMean_rs{s + 1}'] = ak.mean( rseg[0,:,:,2], rseg[0,:,:,0] )
        layer_mean_rs = ak.mean( rseg[0,:,:,3], rseg[0,:,:,0] )
        f_dict[f'xStd_rs{s + 1}'] = math.sqrt(
                ak.mean(
                    (rseg[0,:,:,1] - f_dict[f'xMean_rs{s + 1}'])**2,
                    weight = rseg[0,:,:,0]
                    )
                )
        f_dict[f'yStd_rs{s + 1}'] = math.sqrt(
                ak.mean(
                    (rseg[0,:,:,2] - f_dict[f'yMean_rs{s + 1}'])**2,
                    weight = rseg[0,:,:,0]
                    )
                )
        f_dict[f'layerStd_rs{s + 1}'] = math.sqrt(
                ak.mean(
                    (rseg[0,:,:,3] - layer_mean_rs)**2,
                    weight = rseg[0,:,:,0]
                    )
                )

        # Regions in each segment
        for r, reg in enumerate(rseg):

            # Electron category of region
            cat = ak.concatenate( (reg[0], reg[1]) )
            f_dict[f'eContNHits_x{r + 1}_rs{s + 1}'] = ak.count( cat ) // 4
            if f_dict[f'eContNHits_x{r + 1}_rs{s + 1}'] != 0:
                f_dict[f'eContEnergy_x{r + 1}_rs{s + 1}'] = ak.sum( cat[:,0] )
                f_dict[f'eContXMean_x{r + 1}_rs{s + 1}'] = ak.mean(
                                                                    cat[:,1],
                                                                    cat[:,0]
                                                                    )
                f_dict[f'eContYMean_x{r + 1}_rs{s + 1}'] = ak.mean(
                                                                    cat[:,2],
                                                                    cat[:,0]
                                                                    )
                layer_mean_x_rs = ak.mean( cat[:,3], cat[:,0] )
                f_dict[f'eContXStd_x{r + 1}_rs{s + 1}'] = math.sqrt(
                        ak.mean(
                            (
                                cat[:,1] \
                                - f_dict[f'eContXMean_x{r + 1}_rs{s + 1}']
                                )**2,
                            weight = cat[:,0]
                            )
                        )
                f_dict[f'eContYStd_x{r + 1}_rs{s + 1}'] = math.sqrt(
                        ak.mean(
                            (
                                cat[:,2] \
                                - f_dict[f'eContYMean_x{r + 1}_rs{s + 1}']
                                )**2,
                            weight = cat[:,0]
                            )
                        )
                f_dict[f'eContLayerStd_x{r + 1}_rs{s + 1}'] = math.sqrt(
                        ak.mean(
                            (cat[:,3] - layer_mean_x_rs)**2,
                            weight = cat[:,0]
                            )
                        )

            # Gamma category of region
            cat = ak.concatenate( (reg[0], reg[2]) )
            f_dict[f'gContNHits_x{r + 1}_rs{s + 1}'] = ak.count( cat ) // 4
            if f_dict[f'gContNHits_x{r + 1}_rs{s + 1}'] != 0:
                f_dict[f'gContEnergy_x{r + 1}_rs{s + 1}'] = ak.sum( cat[:,0] )
                f_dict[f'gContXMean_x{r + 1}_rs{s + 1}'] = ak.mean(
                                                                    cat[:,1],
                                                                    cat[:,0]
                                                                    )
                f_dict[f'gContYMean_x{r + 1}_rs{s + 1}'] = ak.mean(
                                                                    cat[:,2],
                                                                    cat[:,0]
                                                                    )
                layer_mean_x_rs = ak.mean( cat[:,3], cat[:,0] )
                f_dict[f'gContXStd_x{r + 1}_rs{s + 1}'] = math.sqrt(
                        ak.mean(
                            (
                                cat[:,1] \
                                - f_dict[f'gContXMean_x{r + 1}_rs{s + 1}']
                                )**2,
                            weight = cat[:,0]
                            )
                        )
                f_dict[f'gContYStd_x{r + 1}_rs{s + 1}'] = math.sqrt(
                        ak.mean(
                            (
                                cat[:,2] \
                                - f_dict[f'gContYMean_x{r + 1}_rs{s + 1}']
                                )**2,
                            weight = cat[:,0]
                            )
                        )
                f_dict[f'gContLayerStd_x{r + 1}_rs{s + 1}'] = math.sqrt(
                        ak.mean(
                            (cat[:,3] - layer_mean_x_rs)**2,
                            weight = cat[:,0]
                            )
                        )

            # Outside category of region
            f_dict[f'oContNHits_x{r + 1}_rs{s + 1}'] = ak.count( cat ) // 4
            if f_dict[f'oContNHits_x{r + 1}_rs{s + 1}'] != 0:
                f_dict[f'oContEnergy_x{r + 1}_rs{s + 1}'] = ak.sum(reg[3,:,0])
                f_dict[f'oContXMean_x{r + 1}_rs{s + 1}'] = ak.mean(
                                                                    reg[3,:,1],
                                                                    reg[3,:,0]
                                                                    )
                f_dict[f'oContYMean_x{r + 1}_rs{s + 1}'] = ak.mean(
                                                                    reg[3,:,2],
                                                                    reg[3,:,0]
                                                                    )
                layer_mean_x_rs = ak.mean( reg[3,:,3], reg[3,:,0] )
                f_dict[f'oContXStd_x{r + 1}_rs{s + 1}'] = math.sqrt(
                        ak.mean(
                            (
                                reg[3,:,1] \
                                - f_dict[f'oContXMean_x{r + 1}_rs{s + 1}']
                                )**2,
                            weight = reg[3,:,0]
                            )
                        )
                f_dict[f'oContYStd_x{r + 1}_rs{s + 1}'] = math.sqrt(
                        ak.mean(
                            (
                                reg[3,:,2] \
                                - f_dict[f'oContYMean_x{r + 1}_rs{s + 1}']
                                )**2,
                            weight = reg[3,:,0]
                            )
                        )
                f_dict[f'oContLayerStd_x{r + 1}_rs{s + 1}'] = math.sqrt(
                        ak.mean(
                            (reg[3,:,3] - layer_mean_x_rs)**2,
                            weight = reg[3,:,0]
                            )
                        )

    # Don't carry around hitinfo after this
    del e_store['ecal_hits']
    del rsegs
''' #Massive irony
