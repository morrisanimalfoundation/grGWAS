Hidden relations between indivisuals is usually uncovered by calculating a genomic relationship matrix (GRM).
GRM is a covariance matrix calculated using the minor allele frequencies (MAFs) of genotyped SNPs assuming a homogeneous population structure.
However, [KING-robust estimation](https://www.ncbi.nlm.nih.gov/pubmed/20926424) is an alternative approach that does not need MAFs. This techqniue is more trusted in mixed-population datasets (unless the parents are from very different populations).

In our study, we can use either techqniues safely. I our tutorial, we will use the KING-robust kinship estimator. 
According to the PLINK2 documentation, the KING kinship coefficients are scaled such that duplicate samples have kinship 0.5. 
Accordingly, first-degree relations (parent-child, full siblings) correspond to ~0.25, second-degree relations correspond to ~0.125, etc. 
It is conventional to use a cutoff of ~0.354 (the geometric mean of 0.5 and 0.25) to screen for monozygotic twins and duplicate samples, ~0.177 to add first-degree relations, etc.
