## Sample-distance and similarity matrices

Assessemnt of similairty between studied individual is useful to identify sample duplications or unexpected outliers. Moreover, it is critical to adjust for hidden population structure.
A pairwise similarity can be calculated by counting the alleles descended from a common ancestor (i.e. Identity by descent; IBD) and/or alleles that are the same, irrespective of whether they are inherited from a recent ancestor (i.e. Identity by state; IBS).

Pedigree information allow us to figure out the expected proportion of shared alleles by descent and thus IBD distances within families. However, the proposed pedigree is not uncommon to be wrong. With genotyping arrays and sequecning, marker information can be used to calculate IBS or IBD similarity matrices and predict the actual pedigree structure 

VanRaden (2008) and Yang et al. (2010) described the genomic relationship matrix (GRM) as a a covariance matrix which uses identical by state (IBS) information scaled by the allele frequencies, as shared rare alleles are more likely to be IBD than common alleles. These methods do not explicitly differentiate between IBD and IBS information [Clark et al. (2012)](https://rune.une.edu.au/web/handle/1959.11/14477).

GRM uses the minor allele frequencies (MAFs) of genotyped SNPs assuming a homogeneous population structure which is not always a true assumption. [KING-robust estimation](https://www.ncbi.nlm.nih.gov/pubmed/20926424) is an alternative approach that does not need MAFs. This techqniue is more trusted in mixed-population datasets (unless the parents are from very different populations).


In the GRLS study, we can use either techqniues safely. Therefore, we will use the KING-robust kinship estimator in this tutorial.
According to the PLINK2 documentation, the KING kinship coefficients are scaled such that duplicate samples have kinship 0.5. 
Accordingly, first-degree relations (parent-child, full siblings) correspond to ~0.25, second-degree relations correspond to ~0.125, etc. 
It is conventional to use a cutoff of ~0.354 (the geometric mean of 0.5 and 0.25) to screen for monozygotic twins and duplicate samples, ~0.177 to add first-degree relations, etc.
