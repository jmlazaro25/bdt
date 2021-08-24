from argparse import ArgumentParser
from os import system

def main():

    # Arguments
    parser = ArgumentParser()
    parser.add_argument(
            'stage',
            default=None,
            help="'train' for training trees or 'eval' for plotting trees"
            )
    parser.add_argument(
            '--bdt',
            default=None,
            help="If stage is eval, specify bdt (without the '_bdt')"
            )
    parser.add_argument(
            '-p',
            dest='procs',
            nargs='+',
            help=('[1:] are assumed to be signal and are hadded together for training'),
            default = ['pn', '0.001', '0.01', '0.1', '1.0']
            )
    args = parser.parse_args()

    ldmx_cmd = 'singularity run --home $PWD $PWD/../ldmx_dev_latest.sif . '

    if args.stage == 'eval':

        for proc in args.procs:
                    
            system(
                    ldmx_cmd + 'hadd ' \
                    + '{bdt}_bdt/evals/{p}_eval.root  \
                       {bdt}_bdt/evals/{p}/*'.format(bdt=args.bdt, p=proc)
                    )

    # Untested - mainframe gives template
    elif args.stage == 'train':

        from glob import glob

        # Hardcoding 'all' path for now because making "all" trees is faster
        # than I thought it would be... and v3/stan because that's all we have
        # NOTE: Might want to move this into batch.py or import it
        train_dir = 'bdt/flats/v3/stan/all/train'

        sfs = []
        for proc in procs[1:]:
            sfs.extend( glob(train_dir + '/' + proc + '/*') )

        system( ldmx_cmd + 'hadd ' \
                + train_dir + '/{}_train.root'.format(procs[0]) + ' ' \
                + train_dir + '/{}'.format(procs[0]) + '/*'
                )
        system( ldmx_cmd + 'hadd ' \
                + train_dir + '/sig_train.root' + ' ' \
                + ' '.join(sfs)
                )

if __name__ == '__main__': main()
