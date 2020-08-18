import Bio
from Bio import SeqIO
from Bio.Seq import Seq

'''
extract reads 50-60 kb
'''

total_records = 0
total_failed = 0
total_passed = 0


with open('nano_chloro.fastq', 'r') as handle:
    with open('50-60k_reads.fastq', 'w') as outfile: 
        for record in SeqIO.parse(handle, 'fastq'):
            total_records += 1
            if len(record.seq) >= 50000 and len(record.seq) <= 60000:
                total_passed += 1
                SeqIO.write(record, outfile, 'fastq')
               
            else:
                total_failed += 1

print(total_records, total_passed, total_failed)
