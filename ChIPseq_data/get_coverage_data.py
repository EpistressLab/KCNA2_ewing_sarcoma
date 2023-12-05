#!/usr/bin/env python

# Import modules
import pysam
import sys
import subprocess as sp
from argparse import ArgumentParser, Action as ArgParseAction

####### Functions #########

def define_options():

        class DebugInfoAction(ArgParseAction):

                def __init__(self, option_strings, dest, **kwargs):
                        super(DebugInfoAction, self).__init__(option_strings, dest, nargs=0, **kwargs)

                def __call__(self, parser, namespace, values, option_string=None):
                        returncode = 0
                        try:
                                get_debug_info()
                        except sp.CalledProcessError as CPE:
                                print("ERROR: {}".format(CPE.output.strip().decode('utf-8')))
                                returncode = CPE.returncode
                        except Exception as e:
                                print("ERROR: {}".format(e))
                                returncode = 1
                        finally:
                                parser.exit(returncode)


        # Argument parsing
        parser = ArgumentParser(description='Create a tsv file with depth coverage of a given genomic region')
        # parser.register('action', 'debuginfo', DebugInfoAction)
        parser.add_argument("-b", "--bam", type=str, required=True,
                help="Path to bam file.")
        parser.add_argument("-c", "--coordinates", type=str, required=True,
                help="Genomic region. Format: chr:start-end. Remember that bam coordinates are 0-based")
        parser.add_argument("-o", "--output", type=str, dest="output", default="coverage_data.tsv",
                help="Path and name of the table output [default=%(default)s]")
        parser.add_argument("--debug-info", action=DebugInfoAction,
                help="Show several system information useful for debugging purposes [default=%(default)s]")
        return parser

def parse_coordinates(c):
        c = c.replace(",", "")
        chr = c.split(":")[0]
        start, end = c.split(":")[1].split("-")
        # Convert to 0-based
        start, end = int(start) - 1, int(end)
        return chr, start, end

def fwd(reads):
    return [read for read in reads if not read.alignment.is_reverse]

def rev(reads):
    return [read for read in reads if read.alignment.is_reverse]

def pp(reads):
    return [read for read in reads if read.alignment.is_proper_pair]

def count_nucleotides(bam_file_path,coordinates):
        chr, start, end = parse_coordinates(coordinates)
        file = pysam.AlignmentFile(bam_file_path, "rb" )
        for col in file.pileup(reference=chr, start=start, end=end):
            chrom = file.getrname(col.tid)
            pos = col.pos
            reads = col.pileups
            yield {'chrom': chrom, 'pos': pos,
                'reads_all': len(reads),
                'reads_fwd': len(fwd(reads)),
                'reads_rev': len(rev(reads)),
                'reads_pp': len(pp(reads)),
                'reads_pp_fwd': len(fwd(pp(reads))),
                'reads_pp_rev': len(rev(pp(reads)))}
        file.close()


if __name__ == "__main__":

        parser = define_options()
        if len(sys.argv)==1:
            parser.print_help()
            sys.exit(1)
        args = parser.parse_args()

        with open(args.output, 'w') as outfile:
            outfile.write("chromosome\tposition\treads_all\treads_fwd\treads_rev\treads_pp\treads_pp_fwd\treads_pp_rev\n")
        
        with open(args.output, 'a') as outfile:
            counts = count_nucleotides(args.bam,args.coordinates)
            for count in counts:
                entry = ""
                for value in count.values():
                    entry += "\t" + str(value)
                entry = entry[1:] + "\n"
                outfile.write(entry)


