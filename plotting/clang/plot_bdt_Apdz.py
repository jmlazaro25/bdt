import argparse
import os
import sys
import time

parser = argparse.ArgumentParser(description = 'Make recoilPt bias plots for a given BDT')

parser.add_argument('--bdtName1',
                    dest = 'bdtName1',
                    default = 'gabrielle',
                    help = 'Name of the BDT version')

parser.add_argument('--bdtName2',
                    dest = 'bdtName2',
                    default = 'segmipv3',
                    help = 'Name of the BDT version that you wish to compare')

parser.add_argument('--bdtName3',
                    dest = 'bdtName3',
                    default = '',
                    help = 'Name of the BDT version that you wish to compare')

parser.add_argument('--flatDir1',
                    dest = 'flatDir1',
                    default = '/nfs/slac/g/ldmx/users/aechavez/ldmx-sw-v2.3.0-w-container/gabrielle_bdt/gabrielle_training/flattrees',
                    help = 'Directory where the flat trees are located')

parser.add_argument('--flatDir2',
                    dest = 'flatDir2',
                    default = '/nfs/slac/g/ldmx/users/aechavez/ldmx-sw-v2.3.0-w-container/segmipv3_bdt/segmipv3_training/flattrees',
                    help = 'Directory where the flat trees of the BDT version that you wish to compare are located')

parser.add_argument('--flatDir3',
                    dest = 'flatDir3',
                    default = '',
                    help = 'Directory where the flat trees of the BDT version that you wish to compare are located')

parser.add_argument('--outDir',
                    dest = 'outDir',
                    default = '/nfs/slac/g/ldmx/users/aechavez/ldmx-sw-v2.3.0-w-container/plots',
                    help = 'Directory where the plots will be saved')

parser.add_argument('--withSelect',
                    dest = 'withSelect',
                    default = 'base',
                    help = 'Selection to make the plots with')

parser.add_argument('--bkgEff',
                    dest = 'bkgEff',
                    default = '0.0001',
                    help = 'Background efficiency corresponding to desired cut')

args = parser.parse_args()

# Fill a list with each process
procs = [0.005, 0.01, 0.05, 'pn']

# Initialize the arguments from the parser
bdtName1 = args.bdtName1
bdtName2 = args.bdtName2
bdtName3 = args.bdtName3
flatDir1 = args.flatDir1
flatDir2 = args.flatDir2
flatDir3 = args.flatDir3
outDir = args.outDir
withSelect = args.withSelect
bkgEff = args.bkgEff

print('Now making plots...')
time.sleep(2)

# Plot recoilPt
if bdtName3 == '' or flatDir3 == '':

    # Pass arguments to the plotting macro
    os.system("root -l -q -b 'plot_bdt_Apdz_0.C+("
        + '"{}", "{}", "{}", "{}", "{}", "{}", {}'.format(bdtName1, bdtName2, flatDir1, flatDir2, outDir, withSelect, bkgEff)
        + ")'"
    )

    # Check whether the plots have saved properly
    for proc in procs:
        if os.path.exists('{}/{}_{}_{}_{}_recoilPt.pdf'.format(outDir, withSelect, bdtName1, bdtName2, proc)):
            pass
        else:
            print('Error: An issue occurred while making the plots!')
            time.sleep(2)
            sys.exit(1)

    print('Plots made successfully!')
    time.sleep(2)

else:

    # Pass the arguments to the plotting macro
    os.system("root -l -q -b 'plot_bdt_recoilPt_1.C+("
        + '"{}", "{}", "{}", "{}", "{}", "{}", "{}", "{}", {}'.format(
            bdtName1,
            bdtName2,
            bdtName3,
            flatDir1,
            flatDir2,
            flatDir3,
            outDir,
            withSelect,
            bkgEff
        )
        + ")'"
    )

    # Check whether the plots have saved properly
    for proc in procs:
        if os.path.exists('{}/{}_{}_{}_{}_{}_recoilPt.pdf'.format(outDir, withSelect, bdtName1, bdtName2, bdtName3, proc)):
            pass
        else:
            print('Error: An issue occurred while making the plots!')
            time.sleep(2)
            sys.exit(1)

    print('Plots made successfully!')
    time.sleep(2)
