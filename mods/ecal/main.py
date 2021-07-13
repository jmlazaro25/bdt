"""
All "Basic" Ecal constants and functions
All space values in mm unless otherwise noted
"""

import numpy as np
from mods import physics

# Gdml values
ecal_front_z = 240.5
sp_thickness = 0.001
clearance = 0.001
ECal_dz = 449.2
ecal_envelope_z = ECal_dz + 1
sp_ecal_front_z = ecal_front_z + (ecal_envelope_z - ECal_dz)/2 - sp_thickness/2 + clearance

sin60 = np.sin(np.radians(60))
module_radius = 85. # Dist form center to midpoint of side
module_side = module_radius/sin60
module_gap = 1.5 # Space between sides of side-by-side mods
cell_radius = 5 # Real circle
cellWidth = 8.7

# DetDescr
LAYER_MASK = 0x3F  # space for up to 64 layers
LAYER_SHIFT = 17
MODULE_MASK = 0x1F  # space for up to 32 modules/layer
MODULE_SHIFT = 12
CELL_MASK = 0xFFF  # space for 4096 cells/module (!)
CELL_SHIFT = 0

ecal_layerZs = ecal_front_z + np.array([7.850,   13.300,  26.400,  33.500,  47.950,
                                        56.550,  72.250,  81.350,  97.050,  106.150,
                                        121.850, 130.950, 146.650, 155.750, 171.450,
                                        180.550, 196.250, 205.350, 221.050, 230.150,
                                        245.850, 254.950, 270.650, 279.750, 298.950,
                                        311.550, 330.750, 343.350, 362.550, 375.150,
                                        394.350, 406.950, 426.150, 438.750        ])


# Arrays holding 68% containment radius per layer for different bins in momentum/angle
radius68_thetalt10_plt500 = [4.045666158618167, 4.086393662224346, 4.359141107602775, 4.666549994726691, 5.8569181911416015, 6.559716356124256, 8.686967529043072, 10.063482736354674, 13.053528344041274, 14.883496407943747, 18.246694748611368, 19.939799900443724, 22.984795944506224, 25.14745829663406, 28.329169392203216, 29.468032123356345, 34.03271241527079, 35.03747443690781, 38.50748727211848, 39.41576583301171, 42.63622296033334, 45.41123601592071, 48.618139095742876, 48.11801717451056, 53.220539860213655, 58.87753380915155, 66.31550881539764, 72.94685877928593, 85.95506228335348, 89.20607201266672, 93.34370253818409, 96.59471226749734, 100.7323427930147, 103.98335252232795]
radius68_thetalt10_pgt500 = [4.081926458777424, 4.099431732299409, 4.262428482867968, 4.362017581473145, 4.831341579961153, 4.998346041276382, 6.2633736512415705, 6.588371889265881, 8.359969947444522, 9.015085558044309, 11.262722588206483, 12.250305471269183, 15.00547660437276, 16.187264014640103, 19.573764900578503, 20.68072032434797, 24.13797140783321, 25.62942209291236, 29.027596514735617, 30.215039667389316, 33.929540248019585, 36.12911729771914, 39.184563500620946, 42.02062468386282, 46.972125628650204, 47.78214816041894, 55.88428562462974, 59.15520134927332, 63.31816666637158, 66.58908239101515, 70.75204770811342, 74.022963432757, 78.18592874985525, 81.45684447449884]
radius68_theta10to20 = [4.0251896715647115, 4.071661598616328, 4.357690094817289, 4.760224640141712, 6.002480766325418, 6.667318981016246, 8.652513285172342, 9.72379373302137, 12.479492693251478, 14.058548828317289, 17.544872909347912, 19.43616066939176, 23.594162859513734, 25.197329065282954, 29.55995803074302, 31.768946746958296, 35.79247330197688, 37.27810357669942, 41.657281051476545, 42.628141392692626, 47.94208483539388, 49.9289473559796, 54.604030254423975, 53.958762417361655, 53.03339560920388, 57.026277390001425, 62.10810455035879, 66.10098633115634, 71.1828134915137, 75.17569527231124, 80.25752243266861, 84.25040421346615, 89.33223137382352, 93.32511315462106]
radius68_thetagt20 = [4.0754238481177705, 4.193693485630508, 5.14209420056253, 6.114996249971468, 7.7376807326481645, 8.551663213602291, 11.129110612057813, 13.106293737495639, 17.186617323282082, 19.970887612094604, 25.04088272634407, 28.853696411302344, 34.72538105333071, 40.21218694947545, 46.07344239520299, 50.074953583805346, 62.944045771758645, 61.145621459396814, 69.86940198299047, 74.82378572939959, 89.4528387422834, 93.18228303096758, 92.51751129204555, 98.80228884380018, 111.17537347472128, 120.89712563907408, 133.27021026999518, 142.99196243434795, 155.36504706526904, 165.08679922962185, 177.45988386054293, 187.18163602489574, 199.55472065581682, 209.2764728201696]
#contRegions = 5 # 4 radii -> 5 regions
contRegions = 3 # Only using 1-3 so don't bother with 4 and 5

##################################################
# For v12 everyg reconstruction
##################################################
mipSiEnergy = 0.130 # MeV
secondOrderEnergyCorrection = 4000./4010.
layerWeights = [
            1.675, 2.724, 4.398, 6.039, 7.696, 9.077, 9.630, 9.630, 9.630, 9.630, 9.630,
            9.630, 9.630, 9.630, 9.630, 9.630, 9.630, 9.630, 9.630, 9.630, 9.630, 9.630,
            9.630, 13.497, 17.364, 17.364, 17.364, 17.364, 17.364, 17.364, 17.364, 17.364,
            17.364, 8.990
            ]

def recE(siEnergy, layer):

    """ Reconstructed energy from sim energy """
    
    return ( (siEnergy/mipSiEnergy) * layerWeights[layer-1] + siEnergy)*\
            secondOrderEnergyCorrection

##################################################
# Get hitID-related info
##################################################

def layer(hit):

    """ Get layerID from ecal hit """
    
    return (hit.getID() >> LAYER_SHIFT) & LAYER_MASK

def module(hit):

    """ Get moduleID from ecal hit """

    return (hit.getID() >> MODULE_SHIFT) & MODULE_MASK

def cell(hit):

    """ Get cellID from ecal hit"""

    return (hit.getID() >> CELL_SHIFT) & CELL_MASK

##################################################
# Feature calculation aiding functions
##################################################

def layerZofHit(hit):

    """ Get layerZ from hit """

    return ecal_layerZs[layer(hit)]

def layerIntercepts(pos,mom,layerZs=ecal_layerZs):

    """ Tuple of projected (x,y)s at each ECal layer """
    
    return tuple( [ physics.projection(pos,mom,z) for z in layerZs ] )

##################################################
# Direct feature calculating functions
##################################################

def base_ecal(f_dict, args, e_store, lq):

    """ Copy base values from LDMX_Events """

    # Get the branch (which is the same and only for all of these)
    ecalVeto = next( iter( args.values() ) )

    f_dict['nReadoutHits']    = ecalVeto.getNReadoutHits()
    f_dict['summedDet']       = ecalVeto.getSummedDet()
    f_dict['summedTightIso']  = ecalVeto.getSummedTightIso()
    f_dict['maxCellDep']      = ecalVeto.getMaxCellDep()
    f_dict['showerRMS']       = ecalVeto.getShowerRMS()
    f_dict['xStd']            = ecalVeto.getXStd()
    f_dict['yStd']            = ecalVeto.getYStd()
    f_dict['avgLayerHit']     = ecalVeto.getAvgLayerHit()
    f_dict['stdLayerHit']     = ecalVeto.getStdLayerHit()
    f_dict['deepestLayerHit'] = ecalVeto.getDeepestLayerHit() 
    f_dict['ecalBackEnergy']  = ecalVeto.getEcalBackEnergy()


def gabrielle_containment(f_dict, args, e_store, lq):

    """ Energy and number containment quantities from input for gabrielle """

    ecalVeto = next( iter( args.values() ) )

    for r in range(contRegions):
        f_dict['electronContainmentEnergy_x{}'.format(r+1)] \
                = ecalVeto.getElectronContainmentEnergy()[r]
        f_dict['photonContainmentEnergy_x{}'.format(r+1)] \
                = ecalVeto.getPhotonContainmentEnergy()[r]
        f_dict['outsideContainmentEnergy_x{}'.format(r+1)] \
                = ecalVeto.getOutsideContainmentEnergy()[r]
        f_dict['outsideContainmentNHits_x{}'.format(r+1)] \
                = ecalVeto.getOutsideContainmentNHits()[r]
        f_dict['outsideContainmentXStd_x{}'.format(r+1)] \
                = ecalVeto.getOutsideContainmentXStd()[r]
        f_dict['outsideContainmentYStd_x{}'.format(r+1)] \
                = ecalVeto.getOutsideContainmentYStd()[r]

def ecal_init(f_dict, args, e_store, lq):

    """
    Calculate/store electron ecalSP and photon targetSP hits
    And layer intercepts for both.
    Used in segcont and mipTracking (and hence segmips)
    """

    for branch_name, branch in args.items():
        if branch_name[:22] == 'TargetScoringPlaneHits': targetSPHits = branch
        if branch_name[:20] == 'EcalScoringPlaneHits': ecalSPHits = branch

    # Get e position and momentum from EcalSP
    e_ecalHit = physics.electronEcalSPHit(ecalSPHits)
    if e_ecalHit != None:
        e_ecalPos, e_ecalP = e_ecalHit.getPosition(), e_ecalHit.getMomentum()

    # Photon Info from targetSP
    e_targetHit = physics.electronTargetSPHit(targetSPHits)
    if e_targetHit != None:
        g_targPos, g_targP = physics.gammaTargetInfo(e_targetHit)
    else:  # Should about never happen -> division by 0 in g_traj
        print('no e at targ!')
        g_targPos = g_targP = np.zeros(3)

    # Get electron and photon trajectories
    e_store['e_traj'] = e_store['g_traj'] = None

    if e_ecalHit != None:
        e_store['e_traj'] = layerIntercepts(e_ecalPos, e_ecalP)

    if e_targetHit != None:
        e_store['g_traj'] = layerIntercepts(g_targPos, g_targP)

    # Store recoil electron momentum magnitude and angle with z-axis
    if e_ecalHit != None:
        e_store['recoilPMag']  = physics.mag(e_ecalP)
    else: e_store['recoilPMag'] = -1.0

    if e_store['recoilPMag'] > 0:
        e_store['recoilTheta'] = physics.angle(e_ecalP, units='radians')
    else: e_store['recoilTheta'] = -1.0
