`parSEQ: ` **ADVANCED PROBE AND RESCUE SEQUENCING FOR
ENHANCED VARIANT RETRIEVAL FROM DNA POOL**

### Read our technical report on parSEQ (Link)


**This python library accompanies our parSEQ technical report. parSEQ is an experimental platform that was developed to ingest bacterial pools of variants -proteins and other- (where each variant is encoded on a plasmid in the bacteria), and deliver sequence-verified, individually separated variants. The process begins with the distribution of bacteria from collective pools into separate clonal wells. After this, parSEQ employs a two-level barcoding system for enhanced multiplexing, incorporating both well-specific and plate-specific barcodes. Once barcoded, samples are merged for Next-Generation Sequencing (NGS), after which the NGS data is analyzed to link each variant with its respective well. Operating within a 384-well plate setup, and supported by automated procedures and Python-based data analysis, parSEQ constitutes an important tool in modern protein engineering processes. While initially dependent on Illumina sequencing for detailed results, recent advancements in Oxford Nanopore (ONT)’s v14 chemistry and improved basecalling algorithms now position Nanopore sequencing as an effective substitute.**

![parSEQ Process Workflow](docs/images/parseq-process-workflow.png)

The figure above describes the parseq process workflow that we have implemented at Adaptyv Biosystems. Refer to our parSEQ technical report section 2.1 (parSEQ Process Description) for a detailed explanation of the process. Note that you can use parSEQ with a different setup. Refer to our parSEQ technical report section 2.2 (Easy Implementation of parSEQ across Labs of Varying Scales) for an overview of the different options available for setting up your own parSEQ process.

### Analyzing parSEQ Results

This library was developed to analyze parSEQ sequencing results. It can be used to analyze sequencing results from both Illumina and Nanopore sequencers. This sequencing library accepts raw sequencing fastq files as input. Each fastq file contains the sequencing output for one 384-well plate. In the figure above, this would be the output from step 3 (Sequencing).
We intended to make this library available for the sequencing output of both Illumina and Nanopore (including different basecalling options for Nanopore), and such, we purposefully start from one fastq-file per plate.

#### Input format

- In case of Illumina paired-end sequencing, please merge the paired-end reads using any opensource tool such as flash or fastp to provide one fastq file per plate. Fastq files should contain only adapter trimmed sequences.
- In case of Nanopore sequencing, please provide one concatenated fastq file per plate post basecalling. Fastq files should contain only adapter trimmed sequences.

Provided fastq files can be zipped or unzipped.

#### Analysis pipeline



