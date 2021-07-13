import math
import numpy as np
from mods.ecal import main as emain

##################################################
# Miscellaneous functions
##################################################

def pos(hit):

    """ Get np.ndarray of hit position """

    return np.array( ( hit.getXPos(), hit.getYPos(), hit.getZPos() ) )

def projection(pos_init, mom_init, z_final):

    """ Project poimt to z_final """

    x_final = pos_init[0] + mom_init[0]/mom_init[2]*(z_final - pos_init[2])
    y_final = pos_init[1] + mom_init[1]/mom_init[2]*(z_final - pos_init[2])

    return (x_final, y_final)

def angle(vec, units, vec2=[0,0,1]):

    """ Angle between vectors (with z by default) """

    if units=='degrees': return math.acos( dot( unit(vec), unit(vec2) ) )*180.0/math.pi
    elif units=='radians': return math.acos( dot( unit(vec), unit(vec2) ) )
    else: print('\nSpecify valid angle unit ("degrees" or "randians")')

def mag(iterable):

    """ Magnitude of whatever """

    return math.sqrt(sum([x**2 for x in iterable]))

def unit(arrayy):

    """ Return normalized np array """

    return np.array(arrayy)/mag(arrayy)

def dot(i1, i2):

    """ Dot iterables """

    return sum( [i1[i]*i2[i] for i in range( len(i1) )] )

def dist(p1, p2):

    """ Distance detween points """

    return math.sqrt(np.sum( ( np.array(p1) - np.array(p2) )**2 ))

def distPtToLine(h1,p1,p2):

    """
    Distance between a point and the nearest point on a line
    defined by endpoints
    """

    return np.linalg.norm(
            np.cross(
                (np.array(h1)-np.array(p1)),
                (np.array(h1)-np.array(p2))
                )
            ) / np.linalg.norm( np.array(p1) - np.array(p2) )

def distTwoLines(h1,h2,p1,p2):

    """ Minimum distance between lines, each line defined by two points """

    e1  = unit( h1 - h2 )
    e2  = unit( p1 - p2 )
    crs = np.cross(e1,e2) # Vec perp to both lines

    if mag(crs) != 0:
        return abs( np.dot( crs,h1-p1) )

    else: # Lines are parallel; need different method
        return mag( np.cross(e1,h1-p1) )

def rotate(point,ang): # move to math eventually

    """ 2D Rotation """

    ang = np.radians(ang)
    rotM = np.array([[np.cos(ang),-np.sin(ang)],
                    [np.sin(ang), np.cos(ang)]])

    return list(np.dot(rotM,point))

##################################################
# e/gamma SP hit info
##################################################

def electronTargetSPHit(targetSPHits):

    """ Get electron target scoringplane hit """

    targetSPHit = None
    pmax = 0
    for hit in targetSPHits:

        if hit.getPosition()[2] > emain.sp_thickness + emain.sp_thickness + 0.5\
                or hit.getMomentum()[2] <= 0 \
                or hit.getTrackID() != 1 \
                or hit.getPdgID() != 11:
            continue

        if mag(hit.getMomentum()) > pmax:
            targetSPHit = hit
            pmax = mag(targetSPHit.getMomentum())

    return targetSPHit

def electronEcalSPHit(ecalSPHits):

    """ Get electron ecal scoringplane hit """

    eSPHit = None
    pmax = 0
    for hit in ecalSPHits:

        if hit.getPosition()[2] > emain.sp_ecal_front_z + emain.sp_thickness/2\
                or  hit.getMomentum()[2] <= 0 \
                or  hit.getTrackID() != 1 \
                or  hit.getPdgID() != 11:
            continue

        if mag(hit.getMomentum()) > pmax:
            eSPHit = hit
            pmax = mag(eSPHit.getMomentum())

    return eSPHit

def electronSPHits(ecalSPHits, targetSPHits):

    """ Get electron target and ecal SP hits """

    ecalSPHit   = electronEcalSPHit(ecalSPHits)
    targetSPHit = electronTargetSPHit(tartgetSPHits)

    return ecalSPHit, targetSPHit

def gammaTargetInfo(eTargetSPHit):

    """ Return photon position and momentum at target """

    gTarget_pvec = np.array([0,0,4000]) - np.array(eTargetSPHit.getMomentum())

    return eTargetSPHit.getPosition(), gTarget_pvec

def gammaEcalSPHit(ecalSPHits):

    """ Get photon ecal scoringplane hit """

    gSPHit = None
    pmax = 0
    for hit in ecalSPHits:

        if hit.getPosition()[2] > emain.sp_ecal_front_z + emain.sp_thickness/2\
                or hit.getMomentum()[2] <= 0 \
                or not (hit.getPdgID() in [-22,22]):
            continue

        if mag(hit.getMomentum()) > pmax:
            gSPHit = hit
            pmax = mag(gSPHit.getMomentum())

    return gSPHit

def elec_gamma_ecalSPHits(ecalSPHits):

    """ Get electron and photon ecal scoringplane hits """

    eSPHit = electronEcalSPHit(ecalSPHits)
    gSPHit = gammaEcalSPHit(ecalSPHits)

    return eSPHit, gSPHit
