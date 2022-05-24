from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
from pathlib import Path
from glob import glob
import os 
import re
import argparse
import initio_utilities as iu

parser = argparse.ArgumentParser(
            prog='creation_of_files_for_Antiviral_Resistance_and_Subtyping_Reports.py',
            usage='%(prog)s -d/--dir <NGS RUN location>',
            description='It takes as an argument the location of the NGS run and stores\
                        an file with 20PC sequences to use for subtyping and another file\
                        with 20PC and 2PC to use for the antiviral resistance genotyping\
                         (Stanford reports)'
                         )
parser.add_argument('--dir', '-d', required=True)
args = parser.parse_args()

path_to_quasibams = Path.cwd()/ args.dir/ 'quasibams'
path_to_refseq = path_to_quasibams.parents[2].joinpath(
                                            'scripts/K03455.1_HXB2.fasta')
run_id = Path(args.dir).stem
path_to_clean_fasta = path_to_quasibams.joinpath('data_to_clean')
INITIO_batch = path_to_clean_fasta.parts[-4]

''' Transfer the clean fasta sequences to the main directory FASTA_sequences '''
FASTA_sequences_dir = path_to_clean_fasta.parents[2] / 'FASTA_sequences'/ run_id
if FASTA_sequences_dir.exists():# == False:
    pass
else:
    os.makedirs(FASTA_sequences_dir)

sequences_for_antiviral_genotyping = []
sequences_for_subtyping = []

for fasta in path_to_clean_fasta.glob('*_all_genes_aligned.fas'):
    path_to_refseq = path_to_clean_fasta.parents[3]/'scripts/K03455.1_HXB2.fasta'
    tempPolymeraseFile = [SeqIO.read( path_to_refseq,'fasta')] 
    SeqRecords = SeqIO.parse(fasta, 'fasta')
    for record in SeqRecords:
        if '.20PC' in record.id:
            seq_id = record.id.replace('.20PC', '_20PC')
            record.id = seq_id
            record.name = ''
            record.description = ''
            record.seq = record.seq.ungap('-').upper()
            sequences_for_subtyping.append(record)
            sequences_for_antiviral_genotyping.append(record)
            tempPolymeraseFile.append(record)
        if '.2PC' in record.id:
            seq_id = record.id.replace('.2PC', '_2PC')
            record.id = seq_id
            record.name = ''
            record.description = ''
            record.seq = record.seq.ungap('-').upper()
            sequences_for_antiviral_genotyping.append(record)
            tempPolymeraseFile.append(record)
    
    tmp_file_name = str(fasta).replace('_all_genes_aligned.fas', 
                                  '_2-20PC_Polymerase_mutation_map.fasta')
    tmp_file_name = tmp_file_name.replace('data_to_clean','tmps')
    SeqIO.write(tempPolymeraseFile, tmp_file_name , 'fasta')
    mafft_cline = MafftCommandline(input= tmp_file_name, maxiterate = 0)
    stdout, stderr = mafft_cline()
    with open(tmp_file_name, "w") as handle:
        handle.write(stdout)

    "create the template to map the mutations in Polymerase using each region"

    polymerase_coordinates = {'Polymerase': [2085, 5096], 
                                'Protease': [2253, 2549],
                                'RT-RNase': [2550, 4229], 
                                'RT': [2550, 3869], 
                                'RNase': [3870, 4229], 
                                'Integrase': [4230, 5096]}
    aln = AlignIO.read(tmp_file_name, 'fasta')
    all_genes_records = [record for record in aln]
    for gene, pos in polymerase_coordinates.items():   
        trimmed_aln, first_pos = iu.trimming_alignment(aln, pos[0], pos[1])
        remaining_length = len(aln[0]) - (first_pos + len(trimmed_aln[0])) 
        gene_seq = ('-'*first_pos + str(trimmed_aln[0].seq) + '-'*remaining_length)
        gene_record = SeqRecord(seq= Seq(gene_seq), 
                                id = (aln[0].id).replace('complete_genome', gene),
                                description= '')
        all_genes_records.append(gene_record)
        name_file = tmp_file_name.replace('tmps', '') #stored in quasibams with tabular files
        SeqIO.write(all_genes_records, Path(name_file),'fasta')

resistance_file_name = Path(FASTA_sequences_dir).joinpath(
                            f'{run_id}_2-20PC_seqs_for_Resistance_report.fasta')
subtyping_file_name = Path(FASTA_sequences_dir).joinpath(  
                            f'{run_id}_20PC_seqs_for_WG_Subtyping.fasta') 
SeqIO.write(sequences_for_antiviral_genotyping, resistance_file_name,'fasta')
SeqIO.write(sequences_for_subtyping, subtyping_file_name , 'fasta')

temporary_files = str(path_to_clean_fasta).replace('data_to_clean','tmps')
shutil.rmtree(Path(temporary_files))


'''2 AND 20 PC FASTA NOT NEEDED ANYMORE'''
for frequency_fasta in path_to_quasibams.glob('*PC.fas'):
    os.remove(frequency_fasta) # 2-20PC