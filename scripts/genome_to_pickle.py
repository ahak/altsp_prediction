import os
import pickle
import sys

chroms_dir = sys.argv[1] 
output_dir = os.path.dirname(chroms_dir.strip('/'))

chrom_fas = [ x for x in os.listdir(chroms_dir) if x.endswith('.fa') ]

hg38 = {}

for fasta in chrom_fas:
    chr_name = fasta.split('.')[0]
    with open(os.path.join(chroms_dir, fasta)) as f:
        s = ""
        print("Reading {}.".format( f.readline().strip() ) )  # reads the sequence name
        for line in f.readlines():
            s = s + line.strip()
    hg38[chr_name] = s


    
with open(os.path.join(output_dir, "hg38.pickle"), 'wb' ) as handle:
    pickle.dump(hg38, handle, protocol = pickle.HIGHEST_PROTOCOL)

    
