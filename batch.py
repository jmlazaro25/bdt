import os
import logging
import subprocess
from glob import glob
from time import sleep
from argparse import ArgumentParser
from configparser import ConfigParser

def main():

    # Syntax
    parser = ArgumentParser()
    parser.add_argument('action')
    parser.add_argument('config')
    parser.add_argument('-n', dest='batch_size', type=int, default=1)
    args = parser.parse_args()

    # Configure logger
    logging.basicConfig(format='[ submitJobs ][ %(levelname)s ]: %(message)s',
            level=logging.DEBUG)

    # Determine how many jobs to submit at a time
    logging.info('Submitting jobs in batches of %s' % args.batch_size)

    # Command that will be used to submit jobs to the batch system
    nfsdir = '/nfs/slac/g/ldmx/users/jmlazaro'
    batch_command = ('bsub '

                     + '-W 1000 '

                     #+ '-q short '
                     #+ '-W 10 '

                     + '-n 3 '
                     + '-R "select[centos7] span[hosts=1]" '
                     + f'singularity run --home $PWD '
                     + f'{nfsdir}/ldmx_pro_vissig.sif . '
                     + 'python3 mainframe.py '
                     + args.action + ' '
                     + args.config + ' '
                     + '--batch '
                    )

    # Read conf file. No options not in conf file accepted for consistency.
    config = ConfigParser()
    config.read(args.config)

    # Build list of complete commands
    job_commands = []
    for proc in config.get('setup', 'processes').replace(' ','').split(','):

        if args.action == 'trees':

            # Shouldn't use, but just because
            trainOrTest = config.get('trees', 'train_or_test')
            if trainOrTest != '':
                print('\nOnly running over {}ing set!!!'.format(trainOrTest))
            sets = {'test', 'train'} if trainOrTest == '' else {trainOrTest}

            for st in sets:

                infiles = sorted( glob(
                    config.get('trees', '{}_indir'.format(st)) +\
                                        '/{}/*.root'.format(proc)
                    ) )

                outdir = config.get( 'trees', '{}_outdir'.format(st) ) \
                                                                  + '/' + proc
                if not os.path.exists( outdir ):
                    os.makedirs( outdir )

                for infile in infiles:

                    job_commands.append(
                            batch_command
                            + '-i ' + infile + ' '
                            + '-o ' + outdir
                            #+ ' -m 100'
                            )

        if args.action == 'eval':

            infiles = sorted( glob(
                    config.get('eval', 'indir') + '/{}/*.root'.format(proc)
                    ) )

            outdir = config.get( 'eval', 'outdir' ) + '/' + proc
            if not os.path.exists( outdir ):
                os.makedirs( outdir )

            for infile in infiles:

                job_commands.append(
                    batch_command
                    + '-i ' + infile + ' '
                    + '-o ' + outdir
                    )

    if args.action == 'train':
        job_commands.append( batch_command )
        #job_commands.append( batch_command + '-m 400') # Only for testing

    # Also for testing (a job at a time)
    #quit( job_commands[0] )

    # Submit them
    for command in job_commands:

        print(command)
        subprocess.Popen(command, shell=True).wait()

        # If a batch of jobs has been submitted, don't submit another batch
        # until all jobs are running.
        if (job_commands.index(command) + 1)%args.batch_size == 0:

            # Initially, wait 10s before submitting other jobs. This lets the
            # batch system catch up so we get an accurate count of pending
            # jobs.
            sleep(15)

            # Check how many jobs are pending
            cjobs = int(
                    subprocess.check_output('bjobs -p | wc -l', shell=True)
                    )
            print('cjobs: %s' % cjobs)
            """ Enough of this i think
            while cjobs != 0:
                logging.info('%s jobs are still pending' % cjobs)
                sleep(30)
                cjobs = int(
                        subprocess.check_output('bjobs -p | wc -l', shell=True)
                        )
                continue
            """

if __name__ == '__main__': main()
