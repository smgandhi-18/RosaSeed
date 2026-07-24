# Evaluation

This directory contains scripts used to evaluate RosaSeed alignment
accuracy against a reference aligner.

## evaluate_alignments.py

Computes:

- Standard mapping accuracy (±50 bp)
- Strict mapping accuracy (±5 bp)
- Structural accuracy
- Sequence-consistent accuracy (SCA)
- False-positive and false-negative read lists

The script supports both **SAM** and **BAM** input through `pysam`.

### Usage

```bash
python3 evaluate_alignments.py \
    bwa.sam \
    rosaseed.sam
```

or

```bash
python3 evaluate_alignments.py \
    bwa.bam \
    rosaseed.bam
```
