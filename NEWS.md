# kmer 1.1.0

* Added "central", "centroid" and "farthest" options for OTU clustering
* Increased sequence hashing speed
* Edits made to vignette to improve clarity
* Enforced a maximum kmer size of 12 for DNA to avoid overflow

# kmer 1.0.2

* Resolved bug in **kcount** where sequence disambiguation was prevented by 
  sapply failing to simplify NULL values.

* Enabled users to manually control compression of amino acid sequences 
  in **kcount** via new 'compress' argument. Thanks to Thomas Shafee for suggestion.

# kmer 1.0.1

* Fixed bug in **otu** that prevented additional arguments being passed to 
nested functions via "dots".


# kmer 1.0.0

* Sumbitted to CRAN 2018-03-05.
