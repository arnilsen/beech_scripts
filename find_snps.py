'''
This script takes the chloroplast alignment of the Dunedin and Nelson
individuals and creates as list of SNPs, excluding indels. It then takes a user defined number
of variable sites to look at and then returns a dictionary with the distance between the
min and max locations of the SNPs in the alignment. The values for each key are the coordinates
in the alignment.

'''

from Bio import SeqIO
import csv
import argparse
import textwrap

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description = 'num variants in window',\
            epilog = textwrap.dedent('''
            Additional information:

            takes the chloroplast alignment of the Dunedin and Nelson
            individuals and creates as list of SNPs, excluding indels. It then takes a user defined number
            of variable sites to look at and then returns a dictionary with the distance between the
            min and max locations of the SNPs in the alignment. The values for each key are the coordinates
            in the alignment.'''))

parser.add_argument('-i', '--input', type = str, metavar = '', required = True, help = 'alignement file')
parser.add_argument('-o', '--output', type = str, metavar = '', required = True, help = 'name of output file')
parser.add_argument('-n', '--num_variables', type = int, metavar = '', required = True, help = 'number of varients in window')
args = parser.parse_args()


def parse_vars(input_file):
    with open(input_file, 'r') as fastfile:
        records = SeqIO.parse(fastfile, 'fasta')
        seq1 = next(records)
        seq2 = next(records)

        positions = []
        #get positions of SNPs
        for i in range(0, len(seq1.seq), 1):
            if seq1[i] != seq2[i] and seq1[i] != '-' and seq2[i] != '-':
                positions.append(i)
        
        num_vars = args.num_variables      
        positions.sort()
        var_window = {}

        count = 0
        for i in range(0, len(positions)-num_vars+1):
            count += 1
            num_var = (str(max(positions[i:i+num_vars])-min(positions[i:i+num_vars])) + '_{}'.format(count))
            var_window[num_var]={}
            var_window[num_var]['coords'] = str(min(positions[i:i+num_vars])) + ':' + str(max(positions[i:i+num_vars]))
        
        return var_window
        #print(var_window)

    
            
def write_out(outfile):
    with open(args.output, 'w') as outfile:
        fnames = ['alignment_len', 'coords']
        var_sites = parse_vars(args.input)
        csv_writer = csv.DictWriter(outfile, fnames, delimiter='\t')
        csv_writer.writeheader()        
        
        for key, value in var_sites.items():
            row = {'alignment_len': key}
            row.update(value)
            csv_writer.writerow(row)
            
write_out(args.input)            
#parse_vars(args.input)