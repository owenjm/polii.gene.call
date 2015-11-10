# polii.gene.call

Rscript for calculating average PolII occupancy and FDR for RNA Pol II DamID datasets.  Based on original algorithms developed by Tony Southall as published in Southall et al. (2013). Dev Cell, 26(1), 101–12. doi:10.1016/j.devcel.2013.05.020.  (Modifications from the original method are described in detail in the sourcecode.)

The script processes datafiles in gatc.gff format, such as those generated by [damidseq_pipeline](https://owenjm.github.io/damidseq_pipeline).

## Requirements

  1. R
  2. The polii.gene.call.r Rscript (make executable if not already and place in your path)
  3. A GFF-formatted list of genes.  A file for release 6 of the Drosophila genome is provided in the archive; most GFF annotation files should also work.  Place this file in an accessible directory and use the --genes.file commandline switch to access it:

		polii.gene.call.r --genes.file=/path/to/my-genes-anotation.gff

## Usage

    polii.gene.call.r [list of gatc.gff files to process]
  
Each file will be processed separately, with the output being two files:

  1. A .csv table of all genes, together with average occupancy and FDR
  2. A plain text list of all genes below the FDR threshold (default is 0.01; change with --FDR= commandline switch)

For a list of all possible commandline options, use

    polii.gene.call.r --help
