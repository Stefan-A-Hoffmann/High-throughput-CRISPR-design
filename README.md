# High-throughput-CRISPR-design

********************************************************************************************************************************************************

CRISPRdesigner.py is intended for the high-throughput design of DNA components for CRISPR/Cas9 mediated editing of individual codons in coding sequences.
It checks for Cas9 targets, whose core recognition is removed by the intended edit. 
If such a target does not exist, it looks for the target closest to the edit and introduces synonymous changes in the donor DNA to abolish cutting by Cas9 

![CRISPRdesigner flowchart](https://user-images.githubusercontent.com/92052871/137486195-fe4875d2-1408-4c4e-a0d9-fd4e42248b8c.jpg)


Required inputs are:

x Flat files of the target genome in forward ("whole_genome.txt") and reverse direction ("whole_genome_rc.txt").
  Supplied here are the respective genomic sequences derived from the R64 genome assembly of Saccharomyces cerevisiae S288C (https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/).

x List of coding sequences of the target genes ("CDS sequences.txt") in the format:

  "Gene name 1"
  "CDS sequence gene 1" \newline
  ...
  
x List of intended edits ("edits.txt") in the format:

  "Gene name" \tab "nucleotide position of edited codon (first base of codon)" \tab "Codon changed for" \newline
  ...
  
x List of targets for each targeted gene ("targets_GeneName.txt") retrieved from CRISPRdirect (http://crispr.dbcls.jp/).
  A script to automatically retrieve these for all to-be-edited genes is supplied here ("CRISPRdirectPOST.py"), taking "edits.txt" and "CDS sequences.txt" as input.
  
The scipt produces "output_file.txt" in the format:

"line number" \tab "Gene name" \tab "nucleotide position of edited codon (first base of codon)" \tab "Codon changed for" \tab "forward guide oligo sequence" \tab "reverse guide oligo sequence" \tab "donor sequence" \tab "forward oligo sequence donor PCR" \tab "reverse oligo sequence donor PCR" \newline
...

The guide oligos have bases added for Golden Gate cloning into pWS082	sgRNA entry vector from the Ellis Lab Yeast CRISPR toolkit.
