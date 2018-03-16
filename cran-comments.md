#Submission 3
This is a patch release addressing a bug in "otu" that prevented 
additional arguments being passed to nested functions via "dots":

### Test environments
 * local ubuntu 16.04.2 x86_64-pc-linux-gnu; R version 3.4.0 
 * travis-ci ubuntu 14.04.5 x86_64-pc-linux-gnu; R 3.4.2
 * winbuilder devel (2018-03-16 r74345)
 
### R CMD check results
There were no ERRORs or WARNINGs.
There was one NOTE:
* checking CRAN incoming feasibility ... NOTE 
Maintainer: ‘Shaun Wilkinson <shaunpwilkinson@gmail.com>’
This is safe to ignore. 

### Downstream dependencies
This is a new package, there are currently no dependencies.  
 

--------------------------------------------------------------------------------

#Submission 2
This is the second submission on 6 Mar 2018 with two changes made on Swetlana's request:

 * Thanks, please write the title in title case:
Fast K-Mer Counting and Clustering for Biological Sequence Analysis
 * Please add a reference for the methods in the 'Description' field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
with no space after 'doi:', 'arXiv:' and angle brackets for auto-linking.
Please fix and resubmit.
Best,
Swetlana Herbrandt

### Test environments
 * local ubuntu 16.04.2 x86_64-pc-linux-gnu; R version 3.4.0 
 * travis-ci ubuntu 14.04.5 x86_64-pc-linux-gnu; R 3.4.2
 * winbuilder devel (2018-03-05 r74345)

### R CMD check results
There were no ERRORs or WARNINGs.
There was one NOTE:
* checking CRAN incoming feasibility ... NOTE 
Maintainer: ‘Shaun Wilkinson <shaunpwilkinson@gmail.com>’
According to several forum discussions this is safe to ignore. 


### Downstream dependencies
This is a new package, there are currently no dependencies.  
  
--------------------------------------------------------------------------------

#Submission 1
This is the original CRAN submission 5 Mar 2018.

### Test environments
 * local ubuntu 16.04.2 x86_64-pc-linux-gnu; R version 3.4.0 
 * travis-ci ubuntu 14.04.5 x86_64-pc-linux-gnu; R 3.4.2
 * winbuilder devel (2018-03-05 r74345)

### R CMD check results
There were no ERRORs or WARNINGs.
There was one NOTE:
    The Title field should be in title case, current version then in title case:
    ‘Fast K-mer Counting and Clustering for Biological Sequence Analysis’
    ‘Fast K-Mer Counting and Clustering for Biological Sequence Analysis’
The original formatting is preferable, but can be changed if needed.

### Downstream dependencies
This is a new package, there are currently no dependencies.
