import sys
import argparse
import warnings
warnings.filterwarnings("ignore")


class Main(object):
    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Python Package for analysing epigenetic data such as ChIP-seq, CUT&RUN, ATAC-seq and CUT&Tag',
            usage='''epimapper <task> [<args>]

Tasks available for using:

    fastqc                      Quality control for raw reads in fastq files

    bowtie2_alignment           Alignment of fastq files to a reference genome

    remove_duplicates           Removes duplicated reads

    fragment_length             Finds and plots fragment lengths in sam files

    filtering                   Filters and converts files from bam to bed

    spike_in_calibration        Spike in normalizes input files

    peak_calling                Calls for enriched regiones in bed files

    heatmap                     Plots heatmaps of enriched regions

    differential_analysis       Preforms differential analysis on input bed files

''')
        parser.add_argument('task', help='Pipeline task to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.task):
            print('***** ERROR: Unrecognized task *****')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.task)()
        
        
    def fastqc(self):
        from .function_scripts.fastqc import set_parser, run
        parser = set_parser(argparse.ArgumentParser(prog='epimapper fastqc',
            description='Quality control for raw reads in fastq files'))
        run(parser.parse_args(sys.argv[2:]))


    def bowtie2_alignment(self):
        from .function_scripts.bowtie2_alignment import set_parser, run
        parser = set_parser(argparse.ArgumentParser(prog='epimapper bowtie2_alignment',
            description='Aligns reads files to reference genome'))
        run(parser.parse_args(sys.argv[2:]))

    def remove_duplicates(self):
        from .function_scripts.remove_duplicates import set_parser, run
        parser = set_parser(argparse.ArgumentParser(prog='epimapper remove_duplicates',
            description='Removes duplicated reads'))
        run(parser.parse_args(sys.argv[2:]))
        
    def fragment_length(self):
        from .function_scripts.fragment_length import set_parser, run
        parser = set_parser(argparse.ArgumentParser(prog='epimapper fragment_length',
            description='Finds and plots the fragment lengths '))
        run(parser.parse_args(sys.argv[2:]))
        
    def filtering(self):
        from .function_scripts.filtering import set_parser, run
        parser = set_parser(argparse.ArgumentParser(prog='epimapper filtering',
            description='Filters and converts files from bam to bed'))
        run(parser.parse_args(sys.argv[2:]))
        
    def spike_in_calibration(self):
        from .function_scripts.spike_in_calibration import set_parser, run
        parser = set_parser(argparse.ArgumentParser(prog='epimapper spike_in_calibration',
            description='Spike in normalizes input files'))
        run(parser.parse_args(sys.argv[2:]))

    def peak_calling(self):
        from .function_scripts.peak_calling import set_parser, run
        parser = set_parser(argparse.ArgumentParser(prog='epimapper peak_calling',
            description='Calls for enriched regiones in bed files'))
        run(parser.parse_args(sys.argv[2:]))
        
    def heatmap(self):
        from .function_scripts.heatmap import set_parser, run
        parser = set_parser(argparse.ArgumentParser(prog='epimapper heatmap',
            description='Plots heatmaps of enriched regions'))
        run(parser.parse_args(sys.argv[2:]))
        
    def differential_analysis(self):
        from .function_scripts.differential_analysis import set_parser, run
        parser = set_parser(argparse.ArgumentParser(prog='epimapper differential_analysis',
            description='Preforms differential analysis on input files'))
        run(parser.parse_args(sys.argv[2:]))

def main():
    Main()

if __name__ == '__main__':
    Main()
