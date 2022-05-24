from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
from glob import glob
from pathlib import Path
import os 
import re
import argparse
import initio_utilities as iu

parser = argparse.ArgumentParser(
            prog='confirmation_of_sequence_cleaning.py',
            usage='%(prog)s -d/--dir <NGS RUN location>',
            description='It takes as an argument the location of the NGS run and\
                         creates a second report to confirm the sequences (all the\
                         20PC and at least POL from 2PC) are in frame'
                         )
parser.add_argument('--dir', '-d', required=True)
args = parser.parse_args()

quasi_parts = [args.dir, 'quasibams']
path_to_quasibams = Path.cwd().joinpath(*quasi_parts)
path_to_refseq = path_to_quasibams.parents[2].joinpath(
                                            'scripts/K03455.1_HXB2.fasta')
run_id = Path(args.dir).stem

os.chdir(path_to_quasibams)
data_to_clean = Path.cwd().joinpath('data_to_clean')
tmps = Path.cwd().joinpath('tmps')
INITIO_batch = Path.cwd().parts[-3]

for fasta in data_to_clean.glob('*.fas'):
    #fasta = INITiO/INITiO2017-2018/Run13_NGS138/
    # quasibams/data_to_clean/RS19002653.1_all_genes_aligned.fas
    # refseq = Path('../../../scripts/K03455.1_HXB2.fasta')
    frequency_check_Pol = [SeqIO.read(path_to_refseq,'fasta')] 
    SeqRecords = SeqIO.parse(fasta, 'fasta')
    for record in SeqRecords:
        if '.20PC' in record.id or '.2PC' in record.id :
            record.name = ''
            record.description = ''
            record.seq = record.seq.ungap('-').upper()
            frequency_check_Pol.append(record)
    tmp_file_name = str(fasta).replace('all_genes_aligned.fas', 'aligned_to_HXB2.fas')
    tmp_file_name = tmp_file_name.replace('data_to_clean','tmps')
    SeqIO.write(frequency_check_Pol, tmp_file_name , 'fasta')
    mafft_cline = MafftCommandline(input= tmp_file_name, maxiterate = 0)
    stdout, stderr = mafft_cline()
    with open(tmp_file_name, "w") as handle:
        handle.write(stdout)

iu.frameshift_check_sequence_locator(run_id, tmps, data_to_clean, check=False)