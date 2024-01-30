from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import os
from pathlib import Path
from itertools import product

def blastn_primer (qry, sub, name, primer, id, cv):
	blastn = NcbiblastnCommandline(query=qry, subject=sub, perc_identity = id, qcov_hsp_perc = cv, task='blastn-short',out=f"/home/sagawa/Documentos/dnaA_project/mit_verify/output2/{primer}/{name}.xml", outfmt=6)
	#blastn = NcbiblastnCommandline(query=qry, subject=sub, task='blastn-short',out=f"/home/sagawa/Documentos/dnaA_project/biopy_blast/output/dnaA_test_foward/{name}.xml", outfmt=6)
	blastn()

def blast_main():
    identity = 94
    cover = 94

    sub_dir = '/home/sagawa/Documentos/dnaA_project/mit_verify/input'
    primer_dir = '/home/sagawa/Documentos/dnaA_project/biopy_blast/loop_dnaA/input/'

    list_sub = os.listdir(sub_dir)
    list_primer = os.listdir(primer_dir)

    for sub in list_sub:
        sub_file = os.path.join(sub_dir, sub)
        out_name = sub.replace('.fasta', '')
        out_path = Path(f'/home/sagawa/Documentos/dnaA_project/mit_verify/output2/{out_name}')
        if not out_path.exists():
            out_path.mkdir()   
    
        for primer in list_primer:
            primer_name = primer.replace('.fasta', '')
            primer_file = os.path.join(primer_dir, primer)
            blastn_primer(primer_file, sub_file, primer_name, out_name, identity, cover)

blast_main()
