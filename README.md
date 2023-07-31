# GWAS of GRLS data

## 1. Setup
## 1.1. Install conda
We are using [conda](https://conda.io/projects/conda/en/stable/index.html) to install softaware needed to run this tutorial. Here are the steps we used to install conda on a **64-bit** computer with a **Linux** system using **Miniconda**. For other operating systems, you can find detailed instructions [here](https://conda.io/projects/conda/en/stable/user-guide/install/index.html)

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

-  Follow the prompts on the installer screens and accept the defaults.
-  Restart the terminal

## 1.2. Create new environment and install softaware packages

```
conda create -n grGWAS
conda activate grGWAS
conda install -c bioconda plink
conda install -c bioconda plink2
conda install -c bioconda bcftools
conda install -c bioconda gcta
conda install -c conda-forge r-base=3.6.3
```
<br>

## 2. Input files
## 2.1. Genotyping data
1.  Affymetrix (thermofisher) Axiom Canine HD Array sets A and B were used to genotype the GRLS dogs. In this tutorial, we have the GRLS genotyping data of each array in a binary PLINK file format. You may use the [PLINK documentation](https://www.cog-genomics.org/plink/1.9/input#bed) to read more about this file format.
    - The files of array A have the prefix "output/setA/export_plink/AxiomGT1.bin"
    - The files of array B have the prefix "output/setB/export_plink/AxiomGT1.bin".
    ---
    **_Note:_** These PLINK files have the gender as predicted by the Axiom genotyping analysis tools.

    ---

2.  A Text file that maps between sample IDs in the genotyping files, biological sample IDs, and the public IDs used in phenotype data files. 
3.  A Text file that has gender information of the dogs as reported by their owners
4.  A test file that has a list of dupliacte samples in the genotyping data 
 
## 2.2. Phenotype data
Morris Animal Foundation [Data Commons](https://datacommons.morrisanimalfoundation.org/) provides open access to most of the data collected by the Golden Retriever Lifetime Study. An overview description of the study data can be found [here](https://datacommons.morrisanimalfoundation.org/node/221). To download data tables, you need to [register](https://datacommons.morrisanimalfoundation.org/user/login?destination=/node/1) at the Data Commons.

In this tutorial, we are using the [Conditions - Neoplasia](https://datacommons.morrisanimalfoundation.org/artisanal_dataset/71) dataset as an example.

<br>


## 3. QC and pre-processing of genotyping data
## 3.1. replicate SNPs  
Both arrays has a number of replicate SNPs which are usefull for QC but also require special attention for proper merging of the files of both arrays. We are using a QC merging mode of PLINK to identify mismatching nonmissing calls between the 2 arrays

```
plink --bfile output/setB/export_plink/AxiomGT1.bin --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --bmerge output/setA/export_plink/AxiomGT1.bin --merge-mode 7 \
      --output-chr 'chrM' --out AxiomGT1.mismatching7
```

Here is the output of this command


> 3339 samples loaded from output/setB/export_plink/AxiomGT1.bin.fam. <br> 
3355 samples to be merged from output/setA/export_plink/AxiomGT1.bin.fam. <br> 
Of these, 17 are new, while 3338 are present in the base dataset. <br> 
547232 markers loaded from output/setB/export_plink/AxiomGT1.bin.bim. <br> 
397685 markers to be merged from output/setA/export_plink/AxiomGT1.bin.bim. <br> 
Of these, 368215 are new, while 29470 are present in the base dataset. <br> 
Warning: Variants 'Affx-206088448' and 'Affx-205939940' have the same position. <br> 
Warning: Variants 'Affx-205859745' and 'Affx-205344060' have the same position. <br> 
Warning: Variants 'Affx-206706550' and 'Affx-205359377' have the same position. <br> 
1460 more same-position warnings: see log file. <br> 
Performing 1-pass diff (mode 7), writing results to AxiomGT1.mismatching7.diff <br> 
98370860 overlapping calls, 97897061 nonmissing in both filesets. <br> 
97593330 concordant, for a concordance rate of 0.996897. <br> 


Let us have a closer look on the mismatching genotypes to see if they are uniformly distributed among the samples or not
```
tail -n+2 AxiomGT1.mismatching7.diff | awk '{print $2,$3}' | sort | uniq -c | sort -k1,1nr > AxiomGT1.mismatching7.diff.samples
```

It seems that two samples show exceptional higher rate of mismatching (S007258 and S019740 have 11676 and 11440 mistmatches respectively. The latter sample is the only sample that was predicted to be female on Array A and male on array B. This likely indicate that these samples had something wrong. Swapping them on array B did not fix the issue)

The mismatching analysis reveals 1463 markrs that has the same position on bith arrays but with different SNP IDs. Futher digging in the array annotation showed that most of these markers are idententical with minor differences in the flanking sequences. To avoid genotyping errors, we will exclude these SNPs from array B.

```
grep "Warning: Variants .* have the same position" AxiomGT1.mismatching7.log | awk -F"'" 'BEGIN{OFS="\n";}{print $2,$4}' > same_pos.arrB.lst
plink --bfile output/setB/export_plink/AxiomGT1.bin --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --exclude same_pos.arrB.lst \
      --make-bed --output-chr 'chrM' --out output/setB/export_plink/AxiomGT1.bin_noSamePos
```
## 3.2. Merging of Array sets A and B Genotyping data
Now, we will merge the genotyping data of both arrays. For SNPs shared between the arrays, the genotypes of array A will overwrite the nonmissing calls in array B. Also, we will exclude the 2 samples with higher rates of non-concordance

```
plink --bfile output/setB/export_plink/AxiomGT1.bin_noSamePos --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --bmerge output/setA/export_plink/AxiomGT1.bin --merge-mode 3 \
      --remove swap_samples.lst \
      --make-bed --output-chr 'chrM' --out AxiomGT1v2.comp_merge
```

The output the merging command gives some useful stastics

> ... <br>
Total genotyping rate in remaining samples is 0.994631. <br>
913984 variants and 3354 samples pass filters and QC. <br>
...


## 3.3. identification and removal of duplicate samples

