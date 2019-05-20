# kmer version 1.1.2

This version provides a bug fix for the kcount function, which wasn't deconstructing
alignments prior to running the kmer counting operation

## Test environments

 * local ubuntu 16.04.2 x86_64-pc-linux-gnu; R version 3.6.0 (2019-04-26)
 * travis-ci ubuntu 14.04.5 x86_64-pc-linux-gnu; R 3.5.2 (2017-01-27)
 * winbuilder R devel  (2019-05-19 r76539)

## R CMD check results

There were no ERRORs WARNINGs or NOTEs.

## Downstream dependencies

- aphid:  OK
- insect: OK


