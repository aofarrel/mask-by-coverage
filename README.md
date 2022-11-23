# mask by coverage

A simple WDL workflow to generate a mask file based on BED coverage via bedtools genomecov.

Input: an array of bam files + integer representation of your minimum coverage
Output: genome coverage for all positions in BEDGRAPH format