# mask by coverage

A simple WDL workflow to generate a mask file based on BED coverage via [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html).

Input: an array of bam files + integer representation of your minimum coverage  
Output: genome coverage for all positions [in BEDGRAPH format](https://genome.ucsc.edu/goldenPath/help/bedgraph.html)  