# Mda5_dermatomyositis_manuscript
Welcome to our Github! Here are the scripts and codes used for our Mda5 dermatomyositis manuscript. Here you'll find the below custom code and installation requirements for:

1.AtoI (AtoG) editing in strand-specific mode by SMARTer seq:
   
   1a. Environment: Python 3
   
   1b. Required software packages:  argparse (v1.1), pandas(v2.3.3), numpy(v2.4.2)
   
2.STAR mapping on HMS O2, with split-read mode, sorting the result bam files by coordinates, and creating the index file

   2a. Environment: HMS O2
   
   2b. Required software packages: gcc(v14.2.0), STAR (2.7.11b), samtools(v1.21)

3. Count TE by class: this utilizes mapped bam files, filter to primary alignments, and count TEs that are differentially expressed (|Log2FC|>1) between patients and HC, and then sum it up by class.
   
   3a. Environment: HMS O2
   
   3b: Required software packages: gcc(v14.2.0), samtools(v1.21), bedtools(2.31.0), python(3.13.1)
   
   3c. Repbase, downloaded directly as a .bed file

   
