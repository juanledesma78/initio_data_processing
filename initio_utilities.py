from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
#from Bio.Data import IUPACData
#from itertools import product
from pathlib import Path
import os
from glob import glob
#import re


def trimming_alignment(BioAln, start_trimming_pos, end_trimming_pos):
    """It takes a Biopython alignment and returns a new alignment trimmed to
    the given starting and ending positions. Gaps (-) introduced in the refseq 
    to get it aligned to the query sequence (insertions) can compromise the 
    final region to trimm. Therefore they are not considered in the numbering
    when trimming is done (if str(BioAln[0].seq[pos]) != '-')
    The positions to trim can be given from a dictionary with several regions to
    trim from the same alignment.
    gene_dict = {'Gag': [790, 2292], 'Pol': [2085, 5096]}
    """

    initial_pos_ungapped = 0
    final_pos_ungapped = 0
    first_pos = 0
    last_pos = 0
    for pos in range(len(BioAln[0].seq)): 
        if initial_pos_ungapped <= (start_trimming_pos-1): # -1 python correction 
            if str(BioAln[0].seq[pos]) != '-':
                first_pos = pos
                initial_pos_ungapped += 1
        if final_pos_ungapped <= end_trimming_pos:
            if str(BioAln[0].seq[pos]) != '-':
                last_pos = pos
                final_pos_ungapped += 1
    aln = BioAln[:, first_pos : last_pos]
    return aln , first_pos


def frameshift_check_sequence_locator(run_id, tmps, data_to_clean, check= True ): 
    """ It takes as arguments RUN id and location of temporary and fasta files.
    It reads the alignments provided in the temporary directory, translates 
    the nt sequences at three diferent frames, aligns their amino acids to the 
    reference sequence and finds out wheter the first one (expected 
    to be in frame after trimming) has the best score. 
    It returns a dictionary with the id of the sequences, the first frame and a comment
    if there is a frameshihf, all of them used to generate a report.
    Thi is an alternative to HIV SEQUENCE LOCATOR 
    https://www.hiv.lanl.gov/content/sequence/LOCATE/locate.html"""

    if check == True:
        seqloc_name = f'FRAMESHIFT_initial_check_{run_id}.txt'
        # CHECK frameshift 
    else:
        seqloc_name = f'FRAMESHIFT_post_cleaning_{run_id}.txt'    

    seqloc_file = open(data_to_clean.joinpath(seqloc_name), 'w')
    seqlocInfo =''

    for aln_file in tmps.glob('*_aligned_to_HXB2.fas'):
        #print(aln_file) #tmps/RS19002656.1_aligned_to_HXB2.fasta
        #sample_name = aln_file.split('/')[-1].replace('_aligned_to_HXB2.fas','')
        sample_name = str(aln_file.name).replace('_aligned_to_HXB2.fas','')
        seqlocInfo +=  f'''
------------------------------------------------------------------------------
---------------- {sample_name} ------------------------------------------------
------------------------------------------------------------------------------

'''

        # (protease: 2253-2549, RT: 2550-3869, RNase: 3870-4229, Integrase:4230-5096)
        # gp120: 6225 → 7757, gp41: 7758 → 8795

        gene_dict = {'Gag': [790, 2292], 'Pol': [2085, 5096], 'Vif': [5041, 5619],
                    'Vpr': [5559, 5850], 'Tat1': [5831, 6045], 'Rev1': [5970, 6045],
                    'Vpu': [6062, 6310], 'Gp160': [6225, 8795], 'Tat2': [8379, 8469],
                    'Rev2': [8379, 8653],'Nef': [8797, 9388]}

        aln = AlignIO.read(aln_file, 'fasta')
        all_genes_records = [record for record in aln] #refseq and 20PC and/or 2PC

        for gene, pos in gene_dict.items():   
            trimmed_aln, first_pos = trimming_alignment(aln, pos[0], pos[1]) # 790-1, 2292)

            ''' This generates a file to map all the genes to the 2-20PC sequences
            Although the genes look aligned the file is not a proper alignment'''
            remaining_length = len(aln[0]) - (first_pos + len(trimmed_aln[0])) 
            gene_seq = ('-'*first_pos + str(trimmed_aln[0].seq) + '-'*remaining_length)
            gene_record = SeqRecord(seq= Seq(gene_seq), 
                                    id = (aln[0].id).replace('complete_genome', gene),
                                    description= '')
            all_genes_records.append(gene_record)
            name_file = str(aln_file).replace('tmps', 'data_to_clean')
            name_file = name_file.replace('_aligned_to_HXB2.fas', 
                                        '_all_genes_aligned.fas')
            SeqIO.write(all_genes_records, Path(name_file),'fasta')

            aaRefSeq = trimmed_aln[0].seq.ungap().translate() # invariable, it should not have gaps, in-frame...
            ntQueySeq = []
            for n in range(1, len(trimmed_aln)):
                ntQueySeq.append(SeqRecord(seq= trimmed_aln[n][:].seq.ungap(), 
                                            id= trimmed_aln[n].id, 
                                             name= ''))
            dict_frame = {}
            for n in range(len(ntQueySeq)):
                scores = []
                frames = []
                frame_analysis = []
                record = ntQueySeq[n]
                for startpos in range(3):
                    #BiopythonWarning: Partial codon, len(sequence) not a multiple of three.
                    codons = []
                    frame = record.seq[startpos:]
                    for n in range(0,len(frame),3):
                        codon = str(frame[n:n+3])
                        if len(codon) == 3:
                            codons.append(codon)
                    frame = Seq(''.join(codons)).translate() # frame multiple of 3
                    frames.append(str(frame))
                    score = pairwise2.align.localxx(frame, aaRefSeq, score_only=True)
                    scores.append(score)
                if scores[0] < scores[1] or scores[0] < scores[2]:
                    frame_analysis.append('Sequence MAY NOT BE in frame')
                dict_frame[record.id] = [frames[0], ''.join(frame_analysis)]
            #return dict_frame 
            for seqid, frameInfo in dict_frame.items():
                AAseq = ''
                for n in range(len(frameInfo[0])):
                    if n !=0 and n%150 == 0:
                        AAseq += '\n'
                        AAseq += frameInfo[0][n]
                    else:
                        AAseq += frameInfo[0][n]
                seqlocInfo += f'{seqid}\t{gene}\t{frameInfo[1]}\n{AAseq}\n'
    
    seqloc_file.write(seqlocInfo)
    seqloc_file.close()


def fasta_selection_by_strings(fasta):
    """It takes a given fasta file and and returns new fasta file(s) containing 
    the sequences which headers match the strings entered in a prompt"""

    search = input("Enter the string(s) to find in the headers separated ONLY by comma: \n")
    searching_list = search.split(',')
    #fasta_input = os.path.join(os.getcwd(),fasta)
    for n in range(len(searching_list)):
        selection = {} # remove duplicates
        for record in SeqIO.parse(fasta, 'fasta'):
            print(record.id) 
            if searching_list[n] in record.id:
                selection[record.id] = SeqRecord(record.seq, record.id, description='')
        if selection.values != '':
            output_file = fasta.replace('.fas',f'_selection_{searching_list[n]}.fasta')
            SeqIO.write(selection.values(), output_file,'fasta')
        #print(selection)


def fasta_selection_by_id(fasta_input, seq_list):
    """It takes a given fasta file with several sequences and returns a new 
    fasta file containing the sequences matching the ids given in txt file"""

    with open(seq_list,'r') as seqid:
        seqidlist = [line.rstrip() for line in seqid ]
    
    selection = {} # remove duplicates
    for record in SeqIO.parse(fasta_input, 'fasta'): #Bio.SeqIO.FastaIO.FastaIterator object
        if record.id in seqidlist:
            selection[record.id] = SeqRecord(record.seq, record.id, description='')
    output_file = fasta_input.replace('.fas','_selection_IDs.fasta')
    SeqIO.write(selection.values(), output_file,'fasta')


def extracting_FASTA_from_file(fasta_input):
    FASTAs = Path(fasta_input.parent).joinpath(f'FASTAs')
    if FASTAs.is_dir() is False:
        os.mkdir(FASTAs)
    for record in SeqIO.parse(fasta_input, 'fasta'):
        SeqIO.write(record, FASTAs/f'{record.id}.fas', 'fasta')

