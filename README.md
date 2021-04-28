# Cas9 targeted enrichment of mobile elements using nanopore sequencing in human genomes

## Guide RNA cleavage-site analysis

Run `MEI_ONTreads_alignment.py` to perform alignments.  
Accepted mobile element names: L1HS, AluYa5, AluYb8, SVA_E, SVA_F
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

## Guide RNA design pipeline

We have included three bash files for guide RNA design: `bash.alu.sh`, `bash.line1.sh`, and `bash.sva.sh` in the `RNA.design.pipelines` folder. 

We also included the consensus sequences for different categories of mobile elements used for the design in the `lib` folder.

Requires:
```
 jellyfish/2.2.8
```

## L1Hs methylation analysis pipeline

We have included two scripts for L1Hs methylation analysis: `non_ref_pipeline.sh` and `reference_piepline.sh` in the `Methylation.pipelines` folder. 

Requires:
```
 nanopolish https://github.com/jts/nanopolish
 methylartist https://github.com/adamewing/methylartist
```

## Citation

* Torrin L. McDonald*,  Weichen Zhou*,  Christopher Castro,  Camille Mumm,  Jessica A. Switzenberg,  Ryan E. Mills,  Alan P. Boyle,
[Cas9 targeted enrichment of mobile elements using nanopore sequencing](https://www.biorxiv.org/content/10.1101/2021.02.10.430605v1), 
bioRxiv, 2021, `https://doi.org/10.1101/2021.02.10.430605`

For PALMER:
* Weichen Zhou, Sarah B Emery, Diane A Flasch, Yifan Wang, Kenneth Y Kwan, Jeffrey M Kidd, John V Moran, Ryan E Mills,
[Identification and characterization of occult human-specific LINE-1 insertions using long-read sequencing technology](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz1173/5680708), 
Nucleic Acids Research, 2019, gkz1173, `https://doi.org/10.1093/nar/gkz1173`

## Contact:

for Nano-Pal and others: arthurz@umich.edu or https://github.com/WeichenZhou

for cleavage-site anlaysis: castrocp@umich.edu or https://github.com/castrocp

for methylation analysis: crmumm@umich.edu or https://github.com/crmumm