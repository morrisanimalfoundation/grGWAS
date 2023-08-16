# GWAS of GRLS data

The Golden Retriever Lifetime Study (GRLS) has a vast amount of phenotypic data collected from veterinarian records, and dog owners’ questioners over years for 3,000+ dogs. In addition, the study has genotyped ~ 1.1 Million genetic markers using Axiom™ Canine Genotyping Array Sets A and B. This reposatory should provide tutorials for handling of the GRLS phenotypic and genotypic data to produce successful GWAS studies

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
conda install -c conda-forge r-base=4.3.1 r-ggplot2=3.4.2 r-gridextra=2.3
```
<br>

## 2. Input files
## 2.1. Genotyping data
1.  Affymetrix (thermofisher) Axiom Canine HD Array sets A and B were used to genotype the GRLS dogs. In this tutorial, we have the GRLS genotyping data of each array in a binary PLINK file format. You may use the [PLINK documentation](https://www.cog-genomics.org/plink/1.9/input#bed) to read more about this file format. The [genotyping analysis notes]() has detailed information on the bioinformatic pipeline used for genotyping and all notes that should be considered before any further analysis.
    - The files of array A have the prefix "output/setA/export_plink/AxiomGT1.bin"
    - The files of array B have the prefix "output/setB/export_plink/AxiomGT1.bin".
    ---
    **_Note:_** These PLINK files have the gender as predicted by the Axiom genotyping analysis tools. Check the [genotyping analysis notes]() for details
    
    ---

2.  A text file that maps between sample IDs in the genotyping files, biological sample IDs, and the public IDs used in phenotype data files. Moreover, the file has gender information of the dogs as reported by their owners: `map_id_sex.tab`
   
## 2.2. Phenotype data
Morris Animal Foundation [Data Commons](https://datacommons.morrisanimalfoundation.org/) provides open access to most of the data collected by the Golden Retriever Lifetime Study. An overview description of the study data can be found [here](https://datacommons.morrisanimalfoundation.org/node/221). To download data tables, you need to [register](https://datacommons.morrisanimalfoundation.org/user/login?destination=/node/1) at the Data Commons.

In this tutorial, we are using the [Conditions - Neoplasia](https://datacommons.morrisanimalfoundation.org/artisanal_dataset/71) dataset as an example.

<br>


## 3. Preparation of Array sets A and B genotyping data
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


The concordance rate of nonmissing genotypes is pretty good (99.7%). However, we still have to handle two issues: 

1.  **Mismatching genotypes:** Let us have a closer look on the mismatching genotypes to see if they are uniformly distributed among the samples or not
    ```
    tail -n+2 AxiomGT1.mismatching7.diff | awk '{print $2,$3}' | sort | uniq -c | sort -k1,1nr > AxiomGT1.mismatching7.diff.samples
    ```

    It seems that two samples show exceptional higher rate of mismatching (`GRLS S007258` and `GRLS S019740` have 11676 and 11440 mistmatches respectively. The latter sample is the only sample that was predicted to be female on Array A and male on array B according to the genotyping analysis notes). This likely indicates that these samples had something wrong. According to the the genotyping analysis notes, swapping the 2 samples on array B did not fix the issue. Therefore, we will exclude both samples from further analysis during the merging step. 

2.  **Variants having the same position:** The mismatching analysis reveals 1463 markrs that has the same position on bith arrays but with different SNP IDs. Futher digging in the array annotation showed that most of these markers are idententical with minor differences in the flanking sequences. To avoid genotyping errors, we will exclude these SNPs from array B.

    ```
    grep "Warning: Variants .* have the same position" AxiomGT1.mismatching7.log | awk -F"'" 'BEGIN{OFS="\n";}{print $2,$4}' > same_pos.arrB.lst
    plink --bfile output/setB/export_plink/AxiomGT1.bin --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --exclude same_pos.arrB.lst \
      --make-bed --output-chr 'chrM' --out output/setB/export_plink/AxiomGT1.bin_noSamePos
    ```
## 3.2. Merging of Array sets A and B Genotyping data
Now, we will merge the genotyping data of both arrays. For SNPs shared between the arrays, the genotypes of array A will overwrite the nonmissing calls in array B. Also, we will **exclude the 2 samples with higher rates of non-concordance**

```
echo "GRLS S007258|GRLS S019740" | tr '|' '\n' > swap_samples.lst
plink --bfile output/setB/export_plink/AxiomGT1.bin_noSamePos --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --bmerge output/setA/export_plink/AxiomGT1.bin --merge-mode 3 \
      --remove swap_samples.lst \
      --make-bed --output-chr 'chrM' --out AxiomGT1v2.comp_merge
```

The output of the merging command gives some useful stastics

> ... <br>
Warning: 1355 het. haploid genotypes present (see AxiomGT1v2.comp_merge.hh ); <br>
Warning: Nonmissing nonmale Y chromosome genotype(s) present; <br>
... <br>
Total genotyping rate in remaining samples is 0.994631. <br>
913984 variants and 3354 samples pass filters and QC. <br>
...

The first warning is is usually caused by male heterozygous calls in the X chromosome pseudo-autosomal region. Looking at `AxiomGT1v2.comp_merge.hh` shows that all the 1355 heterozygous haploid genotypes belong to sample `GRLS S005865`. According to the genotyping analysis notes, this sample is reported male in metadata as well as by the genotyping algorithm on Array B but the algorithm failed to predict its gender on array A. We will handle this sample later during the check for gender accuracy

The second warning indicate that some nonmale samples have Y chromosome genotypes. Most PLINK commands treat these genotypes as missing. This behaviour is intended for safty however it makes it tricky to identify these markers. Here is some workaround to find these markers

```
mkdir -p temp_nonmale
cat AxiomGT1v2.comp_merge.fam | awk '{if($5=="2")print}' > temp_nonmale/female.lst
plink --bfile AxiomGT1v2.comp_merge --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
      --keep temp_nonmale/female.lst --chr y \
      --recode --output-chr 'chrM' --out temp_nonmale/AxiomGT1v2.comp_merge.fe
      
sed -i 's/chrY/chr_notY/' temp_nonmale/AxiomGT1v2.comp_merge.fe.map

plink --file temp_nonmale/AxiomGT1v2.comp_merge.fe --chr-set 38 no-y --allow-no-sex --allow-extra-chr \
      --freq counts \
      --output-chr 'chrM' --out temp_nonmale/AxiomGT1v2.comp_merge.fe.explore
tail -n+2 temp_nonmale/AxiomGT1v2.comp_merge.fe.explore.frq.counts | awk '{if($5!=0 || $6!=0)print}' | wc -l
```
We were able to identify 37 markers nonmissing Y chromosome genotypes in females. These markers could be in the psudoautosomal region but unfortunately we do not have known boundries to those regions in dogs. As mentioned above, it is safe to ignore these markers in PLINK 

## 3.3. Identification and removal of duplicate samples
Plink2 has an efficient function to calculate the KING-robust knickship estimator and filter duplicate samples as well. Duplicate samples have kinship coefficients ~0.5, first-degree relations (parent-child, full siblings) correspond to ~0.25, second-degree relations correspond to ~0.125, etc. Here, we are using the conventional cufoff ~0.354 (the geometric mean of 0.5 and 0.25) to identify and filter duplicate samples

```
plink2 --bfile AxiomGT1v2.comp_merge --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --king-cutoff 0.354 \
       --out AxiomGT1v2.comp_merge.nodup
```

Now, let us compare the files selected to be removed by Plink2 knickship versus the list of samples planned to run in duplicates. This code will identify all the samples planned to run in duplicates then find those remaining after exlcusion of the samples selected by the Plink2 knickship filter. We are expecting one replicate to remain from each group of replicate samples 
```
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]+=1;next}/^Family_ID/{print}{if(a[$3]>1)print}' map_id_sex.tab map_id_sex.tab > dup_samples.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1 FS $2]=1;next}{if(!a[$1 FS $2])print $0}' AxiomGT1v2.comp_merge.nodup.king.cutoff.out.id dup_samples.tab > dup_samples_remain.tab
tail -n+2 dup_samples_remain.tab | cut -f3 | sort | uniq -c | sort -k1,1nr
```
Hear is the list of the remaining duplicates and their counts
> 2 &nbsp;&nbsp;&nbsp; 094-027376 <br>
  1 &nbsp;&nbsp;&nbsp; 094-002188 <br>
  1 &nbsp;&nbsp;&nbsp; 094-002396 <br>
  1 &nbsp;&nbsp;&nbsp; 094-002995 <br>
  . <br>
  . <br>


It seems that 2 duplicates of sample `094-027376` are still remaining! Let's try to calculate their KING knickship to see how similar they are:
```
grep 094-027376 dup_samples.tab > failed_deDup.lst
plink2 --bfile AxiomGT1v2.comp_merge --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --keep failed_deDup.lst --make-king-table \
       --out AxiomGT1v2.comp_merge.failed_deDup
cat AxiomGT1v2.comp_merge.failed_deDup.kin0
```
Here is the output. 
> #FID1 &nbsp;&nbsp;&nbsp; ID1 &nbsp;&nbsp;&nbsp; FID2 &nbsp;&nbsp;&nbsp; ID2 &nbsp;&nbsp;&nbsp; NSNP &nbsp;&nbsp;&nbsp; HETHET &nbsp;&nbsp;&nbsp; IBS0 &nbsp;&nbsp;&nbsp; KINSHIP <br>
GRLS &nbsp;&nbsp;&nbsp; S027376_2 &nbsp;&nbsp;&nbsp; GRLS &nbsp;&nbsp;&nbsp; S027376_1 &nbsp;&nbsp;&nbsp; 879956 &nbsp;&nbsp;&nbsp; 0.0645396 &nbsp;&nbsp;&nbsp; 0.0414475 &nbsp;&nbsp;&nbsp; -0.0559478 <br>

It is obvious that these two samples are unlrelated. Therefore, we will exclude both of them with the duplicate samples selected by the Plink2 knickship filter. These are 120 samples in total so we will end up having 3234 samples in our output PLINK file
```
plink2 --bfile AxiomGT1v2.comp_merge --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --remove <(cat AxiomGT1v2.comp_merge.nodup.king.cutoff.out.id failed_deDup.lst) \
       --make-bed --output-chr 'chrM' --out AxiomGT1v2.comp_merge.deDup
```


## 3.4. Check for gender accuracy & remove samples with wrong gender identities 
According to the [genotyping analysis notes](), there are two samples `GRLS S019740` and `GRLS S005865` that had discordant computed gender on the two arrays. The former was removed already because of the high rate of mismatching genotypes between the 2 arrays (see **3.1. replicate SNPs**). The latter showed high rate heterozygous haploid genotypes despite being reported as male in metadata (see **3.2. Merging of Array sets**). We will discard this sample from further analysis.

Moreover, there are additional 9 samples with concordant gender on both arrays but different from the metadata. We can reproduce this information by known gender in metadata versus the computed gender based on genotyping data
```
echo "Family_ID Individual_ID computed_sex metadata_sex" | tr ' ' '\t' > gender_conflict.lst
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1 FS $2]=$5;next}{if(a[$1 FS $2] && a[$1 FS $2]!=$5)print $1,$2,$5,a[$1 FS $2]}' map_id_sex.tab AxiomGT1v2.comp_merge.deDup.fam >> gender_conflict.lst
```

Wrong gender could be a mistake in the metadata or an indication for sample swap. Therefore, it is safer to exclude these samples from our genotyping data
```
echo "GRLS S005865" | tr ' ' '\t' >> gender_conflict.lst
plink2 --bfile AxiomGT1v2.comp_merge.deDup --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --remove gender_conflict.lst \
       --make-bed --output-chr 'chrM' --out AxiomGT1v2.comp_merge.deDup.sexConfirm
```
The output of our last PLINK command indeicate that our final file has 3224 samples (1615 females, 1609 males; 3224 founders)


## 3.5. Update sample IDs in the genotyping files to match the phenotyping files
The data tables at Morris Animal Foundation [Data Commons](https://datacommons.morrisanimalfoundation.org/) are using grls_ids (e.g. 094-000019) or public_ids (e.g. grlsH764T844). In this step, we will upadate the sample IDs in our genotyping files to match the grls_ids using the `map_id_sex.tab` file that maps between different types of IDs
```
tail -n+2 map_id_sex.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$1,$3}' > grls_id.update.lst 
plink2 --bfile AxiomGT1v2.comp_merge.deDup.sexConfirm --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --update-ids grls_id.update.lst \
       --make-bed --output-chr 'chrM' --out AxiomGT1v2.comp_merge.deDup.sexConfirm.grls_ids
```

We can also do the same to upadate the sample IDs in our genotyping files to match the public_ids
```
tail -n+2 map_id_sex.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$1,$4}' > public_ids.update.lst 
plink2 --bfile AxiomGT1v2.comp_merge.deDup.sexConfirm --chr-set 38 no-xy --allow-no-sex --allow-extra-chr \
       --update-ids public_ids.update.lst \
       --make-bed --output-chr 'chrM' --out AxiomGT1v2.comp_merge.deDup.sexConfirm.public_ids
```

## 4. QC and filtration of genotyping data
## 4.1. Identifcation of 1st degree relatives 

We will use 0.177 (the geometric mean of 0.25 and 0.125) as cutoff for the KING knickship  coeffiecient to identify 1st degree relatives.

```
mkdir -p gwas
plink2 --bfile AxiomGT1v2.comp_merge.deDup.sexConfirm.grls_ids --chr-set 38 no-xy --allow-extra-chr \
       --king-cutoff 0.177 \
       --out gwas/AxiomGT1v2.1st_degree_relatives
```

Let us remove those relatives to avoid inflation of false associations in our GWAS

```
plink2 --bfile AxiomGT1v2.comp_merge.deDup.sexConfirm.grls_ids --chr-set 38 no-xy --allow-extra-chr \
       --remove gwas/AxiomGT1v2.1st_degree_relatives.king.cutoff.out.id \
       --make-bed --output-chr 'chrM' --out gwas/AxiomGT1v2.noRelatives
```

## 4.2. Explore

We use Plink commands to assess heterozygosity, missing genotyping, allele frequency, and hardy-weinberg equilibrium (HWE). The output are the files ending in ".het", ".imiss (Per-individual)/.lmiss (per-variant)", ".frq", and ".hwe"  

```
plink --bfile gwas/AxiomGT1v2.noRelatives --chr-set 38 no-xy --allow-extra-chr \
      --het --missing --freq --hardy   \
      --output-chr 'chrM' --out gwas/AxiomGT1v2.noRelatives.explore

# We currentlky have:
# 913984 variants loaded from .bim file.
# 2519 samples (1279 males, 1240 females) loaded from .fam.
```

## 4.2.1. Transform the heterozygosity output files into a histogram to identify the extreme cases. 

```
awk -v size=0.02 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i){if(i==0) print -1*size,size,a[i]/1;else if(i<0) print (i-1)*size,i*size,a[i]/1;else print i*size,(i+1)*size,a[i]/1 }}'  <(tail -n+2 gwas/AxiomGT1v2.noRelatives.explore.het) > gwas/AxiomGT1v2.noRelatives.explore.het.histo 

cat gwas/AxiomGT1v2.noRelatives.explore.het | awk '{if($6>0.3)print}' >  gwas/AxiomGT1v2.noRelatives.explore.het.high ## high homozygosity

cat gwas/AxiomGT1v2.noRelatives.explore.het | awk '{if(NR==1)print}{if($6<-0.3)print}' >  gwas/AxiomGT1v2.noRelatives.explore.het.low ## low homozygosity
```
We can see at least 2 extreme cases. High heterozygosity may indicate sample contamination, however high homozygosity might just indicate higher inbreeding


## 4.2.2. Transform the missingness output files into a histogram to identify the extreme samples and variants. 

```
awk -v size=0.02 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 gwas/AxiomGT1v2.noRelatives.explore.imiss) > gwas/AxiomGT1v2.noRelatives.explore.imiss.histo 

awk -v size=0.02 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 gwas/AxiomGT1v2.noRelatives.explore.lmiss) > gwas/AxiomGT1v2.noRelatives.explore.lmiss.histo 
```
Missingness per individual shows 13 samples with about ~57% genotyping rates. These samples were genotyped on one array but not the other. While missingness per variant shows 36 extreme variants with about 50% genotyping rates

## 4.2.3. Transform the allele frequency output files into a histogram to identify markers with very low MAF

```
awk -v size=0.01 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 gwas/AxiomGT1v2.noRelatives.explore.frq) > gwas/AxiomGT1v2.noRelatives.explore.frq.histo 

awk -v size=0.001 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 gwas/AxiomGT1v2.noRelatives.explore.frq) > gwas/AxiomGT1v2.noRelatives.explore.frq.histo2 
```
There are 331344 markers that have < 1% MAF, 278461 markers have < 0.1% MAF & 314175 markers have < 0.5% MAF


## 4.2.4. Transform the HWE output files into a histogram to identify markers with significant diviation from HWE

```
awk 'BEGIN{OFS="\t";}{ if($9<1e-10)a["1e-10 or less"]++;
                       else if($9<1e-8)a["1e-8:1e-9"]++; else if($9<1e-9)a["1e-9:1e-10"]++; \
                       else if($9<1e-7)a["1e-7:1e-8"]++; else if($9<1e-6)a["1e-6:1e-7"]++; \
                       else if($9<1e-5)a["1e-5:1e-6"]++; else if($9<1e-4)a["1e-4:1e-5"]++; \
                       else if($9<0.001)a["1e-3:1e-4"]++; else if($9<0.01)a["1e-2:1e-3"]++; } \
                 END { for(i in a) print i,a[i] }'  <(tail -n+2 gwas/AxiomGT1v2.noRelatives.explore.hwe) | sort -g > gwas/AxiomGT1v2.noRelatives.explore.het.histo

cat gwas/AxiomGT1v2.noRelatives.explore.hwe | awk '{if(NR==1)print}{if($9<1e-10 )print}' >  gwas/AxiomGT1v2.noRelatives.explore.hwe.low 
```
1000s of SNPs show significant diviation from HWE, however this is an expected behaviour with selective breeding and artificial selection in dog breeds

## 4.3. Filter
```
echo "Oldies 094-040059" > gwas/het_excess.lst
plink --bfile gwas/AxiomGT1v2.noRelatives --chr-set 38 no-xy --allow-extra-chr \
      --geno 0.05 --mind 0.05 --maf 0.01 \
      --remove gwas/het_excess.lst\
      --make-bed --output-chr 'chrM' --out gwas/AxiomGT1v2.noRelatives.filtered
```
Now, we have 582519 variants and 2505 samples that passed filters and QC.


## 4.4. LD analysis 

## 4.4.1 LD estimation
--r2 reports squared inter-variant allele count correlations
'in-phase' adds a column with in-phase allele pairs
'dprime' adds the absolute value of Lewontin's D-prime statistic
```
plink --bfile gwas/AxiomGT1v2.noRelatives.filtered --chr-set 38 no-xy --allow-extra-chr \
      --r2 'in-phase' 'dprime' \
      --output-chr 'chrM' --out gwas/AxiomGT1v2.noRelatives.filtered.stats

awk -v size=0.01 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($8/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 gwas/AxiomGT1v2.noRelatives.filtered.stats.ld) > gwas/AxiomGT1v2.noRelatives.filtered.stats.ld.r2_histo ## 128035 sequential markers have complete linkage
awk -v size=0.01 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($9/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 gwas/AxiomGT1v2.noRelatives.filtered.stats.ld) > gwas/AxiomGT1v2.noRelatives.filtered.stats.ld.dprime_histo
```

## 4.4.2 LD pruning
```
plink --bfile gwas/AxiomGT1v2.noRelatives.filtered --chr-set 38 no-xy --allow-extra-chr \
      --indep-pairwise 100 10 0.2 \
      --output-chr 'chrM' --out gwas/AxiomGT1v2.noRelatives.filtered.LD_lst ## 567586 of 582519 variants removed.

plink --bfile gwas/AxiomGT1v2.noRelatives.filtered --chr-set 38 no-xy --allow-extra-chr \
      --extract gwas/AxiomGT1v2.noRelatives.filtered.LD_lst.prune.in \
      --make-bed --output-chr 'chrM' --out gwas/AxiomGT1v2.noRelatives.filtered.LD_prune ## 14933 variants remaining
```

## 4.4.3 LD re-evaluation
```
plink --bfile gwas/AxiomGT1v2.noRelatives.filtered.LD_prune --chr-set 38 no-xy --allow-extra-chr \
      --freq --r2 'in-phase' 'dprime' \
      --output-chr 'chrM' --out gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.stats

awk -v size=0.01 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($8/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.stats.ld) > gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.stats.ld.r2_histo

awk -v size=0.01 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($9/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.stats.ld) > gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.stats.ld.dprime_histo

awk -v size=0.001 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.stats.frq) > gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.stats.frq.histo2 ## 278461 markers have < 0.1% MAF & 314175 markers have < 0.5% MAF

paste gwas/AxiomGT1v2.noRelatives.filtered.stats.ld.r2_histo gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.stats.ld.r2_histo | less
```

## 4.5. PCA
```
plink2 --bfile gwas/AxiomGT1v2.noRelatives.filtered.LD_prune --chr-set 38 no-xy --allow-extra-chr \
       --autosome --pca \
       --output-chr 'chrM' --out gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.pca
```

I dentify the principle components that explain most of the variance: 
```
Rscript -e 'val <- read.table("gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.pca.eigenval");'\
'val$varPerc <- val$V1/sum(val$V1);'\
'jpeg(file = "Var_PCs.jpg");'\
'plot( x = seq(1:length(val$varPerc)), y = val$varPerc, type = "o",xlab = "Principle Component", ylab = "Variance explained in %");'\
'dev.off();'
```

![](docs/images/Var_PCs.jpg)<!-- -->



Plot the main principle components:
```
Rscript -e 'require(ggplot2);require(gridExtra);'\
'eigenvec <- read.table("gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.pca.eigenvec", header = TRUE, comment.char="");'\
'eigenvec$X.FID <-  as.factor(eigenvec$X.FID);'\
'plot1 <- ggplot(eigenvec, aes(x = PC1, y = PC2, col = X.FID)) + geom_point() + labs(title = "PCA Plot", x = "PC1", y = "PC2");'\
'plot2 <- ggplot(eigenvec, aes(x = PC3, y = PC4, col = X.FID)) + geom_point() + labs(title = "PCA Plot", x = "PC3", y = "PC4");'\
'plot3 <- ggplot(eigenvec, aes(x = PC1, y = PC3, col = X.FID)) + geom_point() + labs(title = "PCA Plot", x = "PC1", y = "PC3");'\
'plot4 <- ggplot(eigenvec, aes(x = PC2, y = PC3, col = X.FID)) + geom_point() + labs(title = "PCA Plot", x = "PC2", y = "PC3");'\
'combined_plot <- grid.arrange(plot1, plot2, plot3, plot4, nrow = 2);'\
'ggsave("pca_plot.png", combined_plot, width = 12, height = 8, dpi = 400);'
```

![](docs/images/pca_plot.png)<!-- -->


## 5. Pre-process the Skin phenotypic data
## 5.1. Identify affected cases until last year
```
cat phenotypes/conditions_skin.csv | awk 'BEGIN{FS=","}\
    NR==1{print $0;next}{if ($5==1) {\
    if (!($1 in max_year) || $3 > max_year[$1]) {
        max_year[$1] = $3;
        line[$1] = $0;
    }}}END{for (id in line) {print line[id];}}' | tr ',' '\t' > gwas/conditions_skin_lastYear.tab
```

## 5.2. Explore the numbers of cases for each phenotype
```
# Get the number of columns in the header
num_columns=$(head -n 1 "gwas/conditions_skin_lastYear.tab" | awk '{print NF}')
# Loop over each column from the third till the last, and print the sum
for ((i=7; i<=$num_columns; i++)); do
    awk -v col="$i" 'NR==1{dis=$col;next}{sum+=$col} END {print dis,sum}' gwas/conditions_skin_lastYear.tab
done > gwas/conditions_skin.no_of_cases
```

## 5.3. Sample selection
## 5.3.1 Select "cases" of a target phenotype
```
target="atopy"
Tcol=$(head -n 1 "gwas/conditions_skin_lastYear.tab" | awk -v pheno="$target" '{for (i=7; i<=NF; i++){if($i==pheno)print i}}')
awk -v Tcol="$Tcol" 'NR==1{print;next}{if($Tcol=="1")print}' gwas/conditions_skin_lastYear.tab > gwas/conditions_${target}_lastYear.tab

tail -n+2 gwas/conditions_${target}_lastYear.tab | cut -d" " -f1 > gwas/${target}_cases.ids
```

## 5.3.2 Check for co-existing conditions
```
for ((i=7; i<=$num_columns; i++)); do
    awk -v col="$i" 'NR==1{dis=$col;next}{sum+=$col} END {if(sum)print dis,sum}' gwas/conditions_${target}_lastYear.tab
done
```
For now, we will not exclude any samples from cases


## 5.3.3. Identify individuals to be excluded from controls
```
grep 'allergic_reaction\|seasonal_allergy\|angioedema\|facial_edema\|.*_dermatitis\|r_o_atopy\|vaccine_reaction' gwas/conditions_skin.no_of_cases > gwas/${target}.to_be_excluded.lst
while read offTarget;do
  echo $offTarget
  offTcol=$(head -n 1 "gwas/conditions_skin_lastYear.tab" | awk -v pheno="$offTarget" '{for (i=7; i<=NF; i++){if($i==pheno)print i}}')
  awk -v Tcol="$Tcol" -v offTcol="$offTcol" 'NR>1{if($Tcol!="1" && $offTcol=="1")print $1}' gwas/conditions_skin_lastYear.tab
done < <(cut -d" " -f1 gwas/${target}.to_be_excluded.lst) | sort | uniq > gwas/${target}.to_be_excluded.ids
```

## 5.4. Remove the excluded samples from the genotyping dataset
```
cat gwas/${target}.to_be_excluded.ids | grep -Fwf - gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.fam > gwas/${target}.to_be_excluded.samples
plink --bfile gwas/AxiomGT1v2.noRelatives.filtered.LD_prune --chr-set 38 no-xy --allow-extra-chr \
      --remove gwas/${target}.to_be_excluded.samples --maf 0.01 \
      --make-bed --output-chr 'chrM' --out gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.${target}
```
There are 14864 variants and 1885 samples that passed filters and QC.


## 5.5. generate the phenotype file
```
awk 'FNR==NR{a[$1]=1;next}{if(!a[$2])print $1,$2,"1";else print $1,$2,"2";}' gwas/${target}_cases.ids gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.${target}.fam > gwas/${target}.pheno
cat gwas/${target}.pheno | cut -d" " -f3 | sort | uniq -c ## we have 240 cases and 1645 controls
```

## 5.6. Assess distribution of cases in PCA
```
plink2 --bfile gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.${target} --chr-set 38 no-xy --allow-extra-chr \
       --autosome --pca \
       --output-chr 'chrM' --out gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.${target}.pca

Rscript -e 'require(ggplot2);require(gridExtra);'\
'eigenvec <- read.table("gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.atopy.pca.eigenvec", header = TRUE, comment.char="");'\
'pheno <- read.table("gwas/atopy.pheno", header = FALSE, comment.char="");'\
'names(pheno) <- c("FID", "IID", "pheno");'\
'pheno$pheno[pheno$pheno == "1"] <- "Unaffected";pheno$pheno[pheno$pheno == "2"] <- "Affected";'\
'ph_eigenvec <-  merge(pheno, eigenvec, by = "IID");'\
'ph_eigenvec$pheno <-  as.factor(ph_eigenvec$pheno);'\
'plot1 <- ggplot(ph_eigenvec, aes(x = PC1, y = PC2, col = pheno)) + geom_point() + labs(title = "PCA Plot", x = "PC1", y = "PC2");'\
'plot2 <- ggplot(ph_eigenvec, aes(x = PC3, y = PC4, col = pheno)) + geom_point() + labs(title = "PCA Plot", x = "PC3", y = "PC4");'\
'plot3 <- ggplot(ph_eigenvec, aes(x = PC1, y = PC3, col = pheno)) + geom_point() + labs(title = "PCA Plot", x = "PC1", y = "PC3");'\
'plot4 <- ggplot(ph_eigenvec, aes(x = PC2, y = PC3, col = pheno)) + geom_point() + labs(title = "PCA Plot", x = "PC2", y = "PC3");'\
'combined_plot <- grid.arrange(plot1, plot2, plot3, plot4, nrow = 2);'\
'ggsave("pca_plot_atopy.png", combined_plot, width = 12, height = 8, dpi = 400);'
```

![](docs/images/pca_plot_atopy.png)<!-- -->

