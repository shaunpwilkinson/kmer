# kmer version 1.0.2

This is a patch release resolving a bug in *kcount* where sequence 
disambiguation was prevented by sapply failing to simplify NULL values

Original CRAN upload failed for windows with following error:

 * Error in re-building vignettes:
   Error: processing vignette 'kmer-vignette.Rmd' failed with diagnostics:
   DLL 'stringi' not found: maybe not installed for this architecture?
   
The inst/doc directory was deleted for resubmission.

## Test environments

 * local ubuntu 16.04.2 x86_64-pc-linux-gnu; R version 3.4.0 
 * travis-ci ubuntu 14.04.5 x86_64-pc-linux-gnu; R 3.4.4
 * winbuilder devel R version 3.5.0 alpha (2018-04-04 r74529)

## R CMD check results

There were no ERRORs WARNINGs or NOTEs.

## Downstream dependencies

- aphid:  OK


