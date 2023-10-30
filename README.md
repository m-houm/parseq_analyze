`parSEQ: ` **ADVANCED PROBE AND RESCUE SEQUENCING FOR
ENHANCED VARIANT RETRIEVAL FROM DNA POOL**

### Read our technical report on parSEQ (Link)

**This python library accompanies our parSEQ technical report. parSEQ is an experimental platform that was developed to ingest bacterial pools of protein variants (where each variant is encoded on a plasmid in the bacteria), and deliver sequence-verified, individually separated variants. The process begins with the distribution of bacteria from collective pools into separate clonal wells. After this, parSEQ employs a two-level barcoding system for enhanced multiplexing, incorporating both well-specific and plate-specific barcodes. Once barcoded, samples are merged for Next-Generation Sequencing (NGS), after which the NGS data is analyzed to link each variant with its respective well. Operating within a 384-well plate setup, and supported by automated procedures and Python-based data analysis, parSEQ constitutes an important tool in modern protein engineering processes. While initially dependent on Illumina sequencing for detailed results, recent advancements in Oxford Nanopore (ONT)’s v14 chemistry and improved basecalling algorithms now position Nanopore sequencing as an effective substitute.**

![parSEQ Process Workflow](docs/images/parseq-process-workflow.png)

The figure above describes the parseq process workflow that we have implemented at Adaptyv Biosystems. Refer to our parSEQ technical report section 2.1 (parSEQ Process Description) for a detailed explanation of the process. Note that you can use parSEQ with a different setup. Refer to our parSEQ technical report section 2.2 (Easy Implementation of parSEQ across Labs of Varying Scales) for an overview of the different options available for setting up your own parSEQ process.
