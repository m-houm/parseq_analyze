`parSEQ:` **PROBE AND RESCUE SEQUENCING FOR ADVANCED VARIANT RETRIEVAL FROM DNA POOL**
---------
### [Read our technical report on parSEQ](https://www.biorxiv.org/content/10.1101/2023.12.12.571337v1)


**parSEQ is an experimental pipeline that was developed to ingest bacterial pools of protein variants (where each variant is encoded on a plasmid in the bacteria), and deliver sequence-verified, individually separated variants. parSEQ aims to maximizes the capture of protein sequence-function data-sets through a sequence- first-screen-later approach. The process begins with the distribution of bacteria from collective pools into separate clonal wells. After this, parSEQ employs a two-level barcoding system for enhanced multiplexing, incorporating both well-specific and plate-specific barcodes. Once barcoded, samples are merged for Next-Generation Sequencing (NGS), after which the NGS data is analyzed to link each variant with its respective well. Operating within a 384-well plate setup, powered by Next-Generation Sequencing, and supported by automated procedures and Python-based data analysis, parSEQ can sample variant pools with high efficiency, low cost, and fast turnaround times. parSEQ works with both Illumina and Oxford Nanopore Sequencing (ONT).**

**```parseq_analyze``` is a python package that accompanies our parSEQ technical report. It is a python based pipeline to analyze parSEQ NGS results. This pipeline integrates multiple open-source packages and modules, with the main goal of identifying the consensus DNA sequence of every well that parSEQ investigates.**

--------
--------

<h2 style="color:green;">Cite Us</h2>

```bibtex
@article{
author = {Moustafa Houmani, Finlay Peterkin, Gerard Antoun, Louis Fischer, Anissa Hammi},
doi = {10.1101/2023.12.12.571337},
keywords = {NGS, High Throughput Biology, Protein Engineering, Cloud Computing, Automation, Sequence-Function Landscapes},
title = {parSEQ: ADVANCED PROBE AND RESCUE SEQUENCING FOR
ENHANCED VARIANT RETRIEVAL FROM DNA POOL},
url = {https://www.biorxiv.org/content/10.1101/2023.12.12.571337v1},
year = {2023}
}
```


--------
--------

### parSEQ Workflow
![parSEQ Process Workflow](images/parseq-process-workflow.png)

The figure above describes the parseq process workflow that we have implemented at Adaptyv Biosystems. Refer to our parSEQ technical report section 2.1 (parSEQ Process Description) for a detailed explanation of the process. Note that you can use parSEQ with a different setup. Refer to our parSEQ technical report section 2.2 (Easy Implementation of parSEQ across Labs of Varying Scales) for an overview of the different options available for setting up your own parSEQ process.

-------------------
-------------------

### Analyzing parSEQ Results

This package was developed to analyze parSEQ sequencing results. It can be used to analyze sequencing results from both Illumina and ONT sequencers. parseq-analyze accepts raw sequencing fastq files as input. Each fastq file should contain the sequencing output for one 384-well plate. In the figure above, this would be the output from step 3 (Sequencing).
As this library can be used for analysis of parSEQ sequencing output using both Illumina and ONT (including different basecalling options for ONT), we purposefully start from one fastq-file per plate. Check input file format below to make your raw sequencing file compatible with parseq_analyze.

-------------------

#### Raw sequencing file format

- In case of Illumina paired-end sequencing, please **merge the paired-end reads** using any opensource tool such as flash or fastp to provide one fastq file per plate. **Fastq files should have the sequencing adapters trimmed** before starting parseq_analyze. For illumina sequencing, this automatically performed by the sequencing platform (e.g. MiSeq).
- In case of Nanopore sequencing, please provide one concatenated fastq file per plate post basecalling. Fastq files should have the **sequencing adapters trimmed** before starting parseq_analyze. This can be automatically performed on the sequencing platform or using Guppy.
Special Consideration for ONT sequencing: We recommend basecalling using SUP, v4.2.0 Super accuracy basecalling model with Dorado v0.2.5B. This yields higher accuracy basecalling that guarantees higher confidence in well consensus sequence recall. It’s important to note that Dorado necessitates GPU acceleration. If you lack access to GPU acceleration, Dorado can run on Apple’s silicon chips (M1/2 series).

Provided fastq files can be zipped or unzipped.

-------------------

#### Analysis pipeline

parseq_analyze library is intended to be used with a notebook similar to the provided notebook template (parseq_run_notebook). The analysis pipeline is divided into 5 big sections as portrayed in the notebook. Every analysis run will automatically create a json file that tracks all important run information.

> Setup Run
- Ceates a folder with the run name (run directory) in the output directory.
- Copies raw fastq files to a folder named raw_data in the run directory.
- Initializes a json file that will be updated throughout the analysis process to capture all important run information.
- Provides a histograms with stats on the number of sequences, and the average, min, max, and std of sequence length.
- Saves the histograms to /run_directory/fastq_length_histograms/raw_fastq.

> Data Filtering and QC
- Filters raw reads in fastq files for length and phred quality score.
- Outputs a filtered fastq file along with .html and .json quality description files per plate to  /run_directory/fastp_output/plate.
- Outputs length histogram along with length stats to /run_directory/fastq_length_histograms/post_fastp for each plate.
- Updates json file with relevant information.

> Read Demultiplexing
- Uses freebarcodes software to assign a barcode to each of the reads in each of the post_fastp fastq files.
    - Freebarcodes based decoding uses barcoding mapping csv
    - Freebarcodes decoding does barcode search with error provided
    - Outputs text file for each plate with barcode assigned to each read in /run_directory/freebarcodes_output
- Separates reads from freebarcodes decoded file into respective well-level fasta files for each plate in /run_directory/demultiplexed/plate
    - Uses constant regions to get the fwd read of any reverse read, such that the fasta files contain fwd reads only.
    - Separation of freebarcodes decoded files into wells runs in parallel for the plates based on the multiprocessing_cores parameter provided.
- Provides plate level csv files for demultiplexed reads and number of reads per well.
- Updates json file with relevant information.

> Alginment & Consensus

- Aligns every demultiplexed fasta file in  /run_directory/demultiplexed/plate and outputs aligned fasta file to /run_directory/aligned/plate.
    - Alignment algorithm options: mafft, muscle.
    - Alignment runs in parallel for the plates based on the multiprocessing_cores parameter provided.
- Calculates consensus sequence for every well of every plate with the provided consensus_threshold and ambiguous base.
    - Provides plate level csv files for raw consensus and stats /run_directory/consensus/raw_consensus.
    - Consensus runs in parallel for the plates based on the multiprocessing_cores parameter provided.
- Provides a run level csv file for raw consensus and stats under /run_directory.
- Note that consensus sequences output at this step contain the well-level barcodes. They can be trimmed using the trimming and reconstruction oprion below.

> Trimming & Reconstruction (Optional).
- Processes consensus csv files to trim 5' and 3' end based on the trim parameters provided. Trimming is optional.
- Processes consensus csv files to add 5' and 3' end sequences based on the reconstruction parameters provided. Reconstruction is optional.
- Provides plate level trimmed & reconstructed consensus and stats under /run_directory/consensus/trimmed_reconstructed_consensus.
- Provides run level csv file for trimmed & reconstructed consensus and stats under /run_directory.

> Results Visualization
- Creates visualizations based on the consensus sequences.
- Visualizations are based on raw consensus results or trimmed & reconstructed consensus results.
- Creates both plate-level and run-level visualizations
- Visualizations are intended to give a holistic view of the parSEQ run performance.

-------------------
### Installation


You can install `parseq_analyze` directly from `PyPI` via:
```
pip install parseq_analyze
```
The correct packages should be automatically installed, but this is not guaranteed to work as you update your packages/dependencies in the future. Please check the requirements file for dependencies.



