# Cas9 targeted enrichment of mobile elements using nanopore sequencing in human genomes

## Guide RNA cleavage-site analysis

Run scp/MEI_ONTreads_alignment.py, accepts mobile element names: L1HS, AluYa5, AluYb8, SVA_E, SVA_F
```
MEI_ONTreads_alignment.py {mobileElementName} {mobileElementSeq.fasta} {long_reads_input_file.fasta} {outfile.txt}
```

Requires: 
```
import sys 
import Bio
```

## Nano-Pal pipeline for Nanopore reads from Flongle/MinION 

Run for L1Hs
```
bash Nano-Pal.LINE.sh
```

Run for AluYb
```
bash Nano-Pal.AluYb.sh
```

Run for AluYa
```
bash Nano-Pal.AluYa.sh
```

Run for SVA_F
```
bash Nano-Pal.SVA_F.sh
```

Run for SVA_E
```
bash Nano-Pal.SVA_E.sh
```

Requires:
```
 samtools/1.3.1  https://github.com/samtools/samtools
 ncbi-blast++/2.10.0  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
 PALMER https://github.com/mills-lab/PALMER
 minimap2 https://github.com/lh3/minimap2

```
##Contact:
```
for Nano-Pal arthurz@umich.edu
for cleavage-site anlaysis castrocp@umich.edu
```
