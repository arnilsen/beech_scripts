import Bio
from Bio import SeqIO
from Bio.Seq import Seq

'''
remove sequences less than 1000 bp
'''

total_records = 0
total_failed = 0
total_passed = 0

failed_len_list = []

with open('nano_chloro.fastq', 'r') as handle:
    with open('1kb_plus_reads.fastq', 'w') as outfile:
        with open('failed_reads.fastq', 'w') as failed_outfile:
 
            for record in SeqIO.parse(handle, 'fastq'):
                total_records += 1
               
                if len(record.seq) >= 1000:
                    total_passed += 1
                    SeqIO.write(record, outfile, 'fastq')
               
                else:
                    SeqIO.write(record, failed_outfile, 'fastq')
                    total_failed += 1
                    failed_len_list.append(len(record.seq))

max_fail = 0
min_fail = 0

for i in failed_len_list:
    if i > max_fail:
        max_fail = i
    elif i < min_fail:
        min_fail = i

print(sum(failed_len_list)/len(failed_len_list))
print()
print(total_records, total_passed, total_failed)
