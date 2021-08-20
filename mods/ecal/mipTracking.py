import math
import numpy as np
from mods import physics
import mods.ecal.main as emain

##################################################
# Feature calculation aiding functions
##################################################

def nearPhotonInf(trackingHitList, g_trajectory):

    """
    First layer with hit within one cell of the photon trajectory AND
    Number of hits within one cell of the photon trajectory
    """

    layer = 33
    n = 0
    for hit in trackingHitList:

        # Near the photn trajectory
        if physics.dist(
                            physics.pos(hit)[:2],
                            g_trajectory[ emain.layer( hit ) ]
                        ) < emain.cellWidth:
            n += 1

            # Earliest layer
            if emain.layer( hit ) < layer:
                layer = emain.layer( hit )

    return layer, n

# Based on previous python; All of v9 analysis done with this
def straightTracks(hitlist, etraj_ends, ptraj_ends, mst=2):

    """ Straight-track mip tracking algorithm """

    strtracklist = []   # Initialize output

    for hit in hitlist:  #Go through all hits, starting at the back of the ecal

        track = [hit]
        currenthit = hit  #Stores "trailing" hit in track being constructed
        possibleNeigh = False

        for h in hitlist:

            if h.getZPos() == currenthit.getZPos():
                possibleNeigh = True  #Optimization
                continue
            if not possibleNeigh:  continue
            if currenthit.getZPos() - h.getZPos() > 25:  #Optimization
                possibleNeigh = False
                continue

            neighFound = (
                    (
                        emain.layer(h) == emain.layer( currenthit ) - 1 \
                        or emain.layer(h) == emain.layer( currenthit ) - 2 \
                    ) \
                    and h.getXPos() == currenthit.getXPos() \
                    and h.getYPos() == currenthit.getYPos()
                )

            if neighFound:
                track.append(h)
                currenthit = h

        # Too short
        if len(track) < mst: continue

        # If it's exactly the min, it has to be very close to ptraj
        if len(track) == mst:
            for hitt in track:
                if physics.distPtToLine( physics.pos(hitt),
                        ptraj_ends[0], ptraj_ends[1] ) > 8:
                    break
                continue

        # Check that the track approaches the photon's and not the electron's
        trk_s = np.array( (
                            track[ 0].getXPos(),
                            track[ 0].getYPos(),
                            track[ 0].getZPos()
                            ) )
        trk_e = np.array( (
                            track[-1].getXPos(),
                            track[-1].getYPos(),
                            track[-1].getZPos()
                            ) )
        closest_e = physics.distTwoLines(
                                            etraj_ends[0], etraj_ends[1],
                                            trk_s, trk_e
                                            )
        closest_p = physics.distTwoLines(
                                            ptraj_ends[0], ptraj_ends[1],
                                            trk_s, trk_e
                                            )

        if closest_p > emain.cellWidth and closest_e < emain.cellWidth:
            continue

        # Remove hits in current track from further consideration
        for h in track:
            hitlist.remove(h)

        # Append track to track list
        strtracklist.append(track)

    # Combine tracks that should be consecutive
    # NOTE: Should maybe do this eariler in case 2 len=2 tracks
    # add up to a passing 4
    strtracklist.sort(key=lambda h: hit.getZPos(), reverse=True)

    currentInd = 0
    while currentInd < len(strtracklist):

        trk = strtracklist[currentInd]
        tmpInd = currentInd+1
        mergeFound = False

        # Search for track compatible with current one
        while tmpInd < len(strtracklist) and not mergeFound:
            trk_ = strtracklist[tmpInd]
            trk_e = np.array( (track[-1].getXPos(), track[-1].getYPos(),
                                                    track[-1].getZPos() ) )
            trk_s = np.array( (track[ 0].getXPos(), track[ 0].getYPos(),
                                                    track[ 0].getZPos() ) )

            # If head+tail are w/in one cell of each other
            if physics.dist( trk_e, trk_s ) < emain.cellWidth:
                for hit in trk_:
                    trk.append(hit) # add hits of one to other
                strtracklist.remove(trk_) # remove one as seperate track
                mergeFound = True
            tmpInd += 1
        if not mergeFound:
            currentInd += 1

    return len(strtracklist)

# Based on C++ Analyzer
def nStraightTracks_c(trackingHitList, e_traj_ends, g_traj_ends):

    """ Straight-track mip tracking algorithm (c++ based version) """

    nTracks = 0

    # Seed a track with each hit
    iHit = 0
    while iHit < len(trackingHitList):
        track = 34*[999]
        track[0] = iHit
        currentHit = iHit
        trackLen = 1

        # Search for hits in next two layers
        jHit = 0
        while jHit < len(trackingHitList):

            if trackingHitList[jHit].layer \
                    == trackingHitList[currentHit].layer \
                or trackingHitList[jHit].layer \
                    > trackingHitList[currentHit].layer + 2:
                jHit += 1 # Don't keep checking this hit over and over again
                continue # Continue if not in the right range

            # If it's also directly behind the current hit, add it to the
            # current track
            if trackingHitList[jHit].pos[:1] \
                    == trackingHitList[currentHit].pos[:1]:

                track[trackLen] = jHit
                currentHit = jHit # Update end of track
                trackLen += 1

            jHit += 1 # Move j along

        # Confirm if track is valid
        if trackLen >= 2: # Set min track length

            # Make sure the track is near the photon trajectory and away from
            # the electron
            closest_e = physics.distTwoLines(
                                    trackingHitList[ track[0] ].pos,
                                    trackingHitList[ track[trackLen-1] ].pos,
                                    e_traj_ends[0], e_traj_ends[1]
                                    )
            closest_g = physics.distTwoLines(
                                    trackingHitList[ track[0] ].pos,
                                    trackingHitList[ track[trackLen-1] ].pos,
                                    g_traj_ends[0], g_traj_ends[1]
                                    )
            if closest_g > physics.cellWidth \
                    and closest_e < 2*physics.cellWidth:
                iHit += 1; continue
            if trackLen < 4 and closest_e > closest_g:
                iHit += 1; continue

            # If valid track is found, remove hits in track from hitList
            for kHit in range(trackLen):
                trackingHitList.pop( track[kHit] - kHit)

            # nStraightTracks++
            nTracks += 1

            # Decrease iHit because the *current" seed will have been removed
            iHit -= 1

        iHit += 1 # Move iHit along

        # Possibley merge tracks later

    return nTracks

##################################################
# Direct feature calculating functions
##################################################

def mipTracking_init(f_dict, args, e_store, lq):

    """ Init prereqs for mipTracking """
    # Some of this could be other init funcs, but for now, whatever

    e_store['trackingHitList'] = []

    if e_store['e_traj'] != None and e_store['g_traj'] != None:

        # Create arrays marking start and end of each trajectory
        e_store['e_traj_ends'] = [
                                    np.array([
                                        e_store['e_traj'][0][0],
                                        e_store['e_traj'][0][1],
                                        emain.ecal_layerZs[0]
                                        ]),
                                    np.array([
                                        e_store['e_traj'][-1][0],
                                        e_store['e_traj'][-1][1],
                                        emain.ecal_layerZs[-1]
                                        ])
                                    ]
        e_store['g_traj_ends'] = [
                                    np.array([
                                        e_store['g_traj'][0][0],
                                        e_store['g_traj'][0][1],
                                        emain.ecal_layerZs[0]
                                        ]),
                                    np.array([
                                        e_store['g_traj'][-1][0],
                                        e_store['g_traj'][-1][1],
                                        emain.ecal_layerZs[-1]])
                                    ]

    else:

        # Electron trajectory is missing so all hits in Ecal are okay to use
        # Pick trajectories so they won't restrict tracking, far outside Ecal

        e_store['e_traj_ends'] = [
                                    np.array([999 ,999 ,0]),
                                    np.array([999 ,999 ,999 ])
                                    ]
        e_store['g_traj_ends'] = [
                                    np.array([1000,1000,0]),
                                    np.array([1000,1000,1000])
                                    ]

    # Territory setup
    e_store['gToe'] = physics.unit(
                        e_store['e_traj_ends'][0] - e_store['g_traj_ends'][0])
    e_store['origin'] = e_store['g_traj_ends'][0] \
                        + 0.5*emain.cellWidth*e_store['gToe']

def epSepAndDot(f_dict, args, e_store, lq):

    """
    Distance between electron and photon at ecal face
    And dot product of normalized electron and photon trajectories
    """

    if e_store['e_traj_ends'][0][0] == 999:

        f_dict['epSep'] = 10.0 + 1.0 # Don't cut on these in this case
        f_dict['epDot'] = 3.0 + 1.0

    else:

        f_dict['epSep'] = physics.dist(
                                        e_store['e_traj_ends'][0],
                                        e_store['g_traj_ends'][0]
                                        )
        f_dict['epDot'] = physics.dot(
                physics.unit(
                    e_store['e_traj_ends'][1] - e_store['e_traj_ends'][0]
                    ),
                physics.unit(
                    e_store['g_traj_ends'][1] - e_store['g_traj_ends'][0]
                    )
                )

def mipTracking_collect(f_dict, args, e_store, ecalRecHit):

    """
    Collect hits for mipTracking
    (outside electron region or electron missing)
    """

    if ecalRecHit.getEnergy() <= 0: return

    # Reused
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

    # If hit's distance to electron greater than RoC at given layer or e-miss
    if distance_e_traj >= e_store['e_radii'][layer] or distance_e_traj == -1.0:
        e_store['trackingHitList'].append(ecalRecHit)

def fullTerritories(f_dict, args, e_store, ecalRecHit):

    """ Full territory selections """

    hitPrime = physics.pos(ecalRecHit) - e_store['origin']

    if np.dot(hitPrime, e_store['gToe']) > 0:
        f_dict['fullElectronTerritoryHits'] += 1
    else:
        f_dict['fullPhotonTerritoryHits'] += 1

def territories(f_dict, args, e_store, lq):

    """ Partial territory and ratio features """

    # Territories limited to trackingHitList

    if e_store['e_traj'] != None:

        for hit in e_store['trackingHitList']:

            hitPrime = physics.pos(hit) - e_store['origin']

            if np.dot(hitPrime, e_store['gToe']) > 0:
                f_dict['electronTerritoryHits'] += 1
            else: f_dict['photonTerritoryHits'] += 1

    else:
        f_dict['photonTerritoryHits'] = f_dict['nReadoutHits']
        f_dict['TerritoryRatio'] = 10
        f_dict['fullTerritoryRatio'] = 10

    if f_dict['electronTerritoryHits'] != 0:

        f_dict['TerritoryRatio'] \
                = f_dict['photonTerritoryHits']/f_dict['electronTerritoryHits']

    if f_dict['fullElectronTerritoryHits'] != 0:

        f_dict['fullTerritoryRatio'] \
                = f_dict['fullPhotonTerritoryHits']\
                    /f_dict['fullElectronTerritoryHits']

def nearPhotonInfo(f_dict, args, e_store, lq):

    """ Just call nearPhotonInf if possible"""

    if e_store['g_traj'] != None:

        f_dict['firstNearPhLayer'], f_dict['nNearPhHits'] \
            = nearPhotonInf( e_store['trackingHitList'], e_store['g_traj'] )

def tracks(f_dict, args, e_store, lq):

    """ Find desired MIP tracks """

    f_dict['straight4'] = straightTracks(
                                                e_store['trackingHitList'],
                                                e_store['e_traj_ends'],
                                                e_store['g_traj_ends'],
                                                mst = 4
                                            )
