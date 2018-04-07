# kmer version 1.0.2

This is a patch release resolving a bug in *kcount* where sequence 
disambiguation was prevented by sapply failing to simplify NULL values

## Test environments

 * local ubuntu 16.04.2 x86_64-pc-linux-gnu; R version 3.4.0 
 * travis-ci ubuntu 14.04.5 x86_64-pc-linux-gnu; R 3.4.2
 * winbuilder devel (2018-03-05 r74345)

## R CMD check results

There were no ERRORs or WARNINGs.
There was one NOTE:
* checking CRAN incoming feasibility ... NOTE 
Maintainer: ‘Shaun Wilkinson <shaunpwilkinson@gmail.com>’
According to forum discussions this is safe to ignore. 


## Downstream dependencies

- aphid:  OK


