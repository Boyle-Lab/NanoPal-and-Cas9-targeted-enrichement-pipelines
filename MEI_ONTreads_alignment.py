#!/usr/bin/env python3

###########################
#  Run this as:
# 'script_name.py {mobileElementName} {mobileElementSeq.fasta} {long_reads_input_file.fasta} {outfile.txt}'
#  Accepted mobile element names: L1HS, AluYa5, AluYb8, SVA_E, SVA_F
##########################

import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align

mobile_element_name = sys.argv[1]
mobile_element_fasta = sys.argv[2]
read_fasta = sys.argv[3]
outfile = sys.argv[4]

########################################################
# Store mobile element sequence and reverse complement #
########################################################
mobile_element = SeqIO.read(mobile_element_fasta, "fasta")
mobile_element_seq = mobile_element.seq
rev_comp_mobile_element_seq = mobile_element_seq.reverse_complement()

############################################
# Store length of mobile element sequences #
############################################

L1HS_len = 6059
AluYa5_len = 311
AluYb8_len = 318
SVA_E_len = 1382
SVA_F_len = 1375

if mobile_element_name == "L1HS":
    mobile_element_len = L1HS_len

elif mobile_element_name == "AluYa5":
    mobile_element_len = AluYa5_len

elif mobile_element_name == "AluYb8":
    mobile_element_len = AluYb8_len

elif mobile_element_name == "SVA_E":
    mobile_element_len = SVA_E_len

elif mobile_element_name == "SVA_F":
    mobile_element_len = SVA_F_len

##################
# Create Aligner #
##################

aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 5
aligner.mismatch_score = -5
aligner.open_gap_score = -10
aligner.extend_gap_score = -5

##############################
# Align reads to L1 sequence #
##############################
with open(outfile, "w") as results_out:
    i = 1 # Keep track of how many reads have been processed
    
    for seq in SeqIO.parse(read_fasta, "fasta"):
        if i % 50 == 0:
            print(str(i) + " reads processed")
        results_out.write(seq.id + "\n")

        # Use Bio.Align package to do pairwise alignments

        forward_alignments = aligner.align(seq.seq[:80], mobile_element_seq)
        if not any(a.score >= 100 for a in forward_alignments):
            results_out.write("No alignments with score >= 100 on the forward strand" + "\n")
        else:
            if ((forward_alignments[0].aligned)[0][0][0]) == 0:
                results_out.write("Forward alignment score: " + "\t" + str(forward_alignments[0].score) + "\n")
                results_out.write(str(forward_alignments[0].aligned) + "\n")
                results_out.write("Forward alignment starting mobile element base:" + "\t" + str((forward_alignments[0].aligned)[1][0][0] + 1) + "\t" + mobile_element_name + "\n")
            else:
                results_out.write("Forward strand alignment does not begin at start of read" + "\n")

        rev_comp_alignments = aligner.align(seq.seq[:80], rev_comp_mobile_element_seq)
        if not any(a.score >= 100 for a in rev_comp_alignments):
            results_out.write("No alignments with score >= 100 on the reverse strand" + "\n")
            results_out.write("\n")
        else:
            if ((rev_comp_alignments[0].aligned)[0][0][0]) == 0:
                results_out.write("Reverse strand alignment score: " + "\t" + str(rev_comp_alignments[0].score) + "\n")
                results_out.write(str(rev_comp_alignments[0].aligned) + "\n")
                results_out.write("Reverse alignment starting mobile element base:" + "\t" + str( mobile_element_len - ((rev_comp_alignments[0].aligned)[1][0][0])) + "\t" + mobile_element_name + "\n")
                results_out.write("\n")  #subtracting the first alignment base from the size of the element to get the starting base with respect to the forward mobile element sequence
            else:
                results_out.write("Reverse strand alignment does not begin at start of read" + "\n")
                results_out.write("\n")
                
        i += 1
  
