## Admixture Analysis
## 1. install
```
conda install -c bioconda admixture=1.3.0
```

## 2. remove non-autosomal markers and change chr names to intigers
```
plink --bfile gwas/AxiomGT1v2.noRelatives.filtered.LD_prune --chr-set 38 no-xy --allow-extra-chr \
      --autosome \
      --make-bed -out gwas/AxiomGT1v2.noRelatives.filtered.LD_prune.autosomes
```

## 3. Run Admixture in the cross validation mode to find out the best number of ancestries
```
cd gwas
for k in {1..10};do echo $k;
  admixture -j8 --cv AxiomGT1v2.noRelatives.filtered.LD_prune.autosomes.bed $k | tee log${k}.out
done
grep "^CV error" log*.out | cut -d":" -f1,3 | sed 's/log//;s/\.out://' | sort -k1,1n > admix_cv.txt

Rscript -e 'val <- read.table("admix_cv.txt");'\
'jpeg(file = "admix_cv.jpg");'\
'plot( x = val$V1, y = val$V2, type = "o",xlab = "K Ancestry", ylab = "CV error");'\
'dev.off();'
```

![](images/admix_cv.jpg)<!-- -->


There is an output file for each parameter set: Q (the ancestry fractions), and P (the allele frequencies of the inferred ancestral populations). 

## 4. plot the Q estimates
```
sort -k1,3nr AxiomGT1v2.noRelatives.filtered.LD_prune.autosomes.3.Q > AxiomGT1v2.noRelatives.filtered.LD_prune.autosomes.3.Q.sorted
Rscript -e 'tbl=read.table("AxiomGT1v2.noRelatives.filtered.LD_prune.autosomes.3.Q.sorted");'\
'jpeg(file = "ancestries.jpg");'\
'barplot(t(as.matrix(tbl)), col=rainbow(3),xlab="Individual #", ylab="Ancestry", border=NA);'\
'dev.off();'
```

![](images/ancestries.jpg)<!-- -->


## 5. plot the Q estimates of a selected subset of individuals
```
tail -n+2 AxiomGT1v2.noRelatives.filtered.LD_prune.pca.eigenvec | awk '{if($3<0.02 && $4>0.03)print $1,$2}' > pc2_extreme
paste AxiomGT1v2.noRelatives.filtered.LD_prune.autosomes.3.Q  AxiomGT1v2.noRelatives.filtered.LD_prune.autosomes.fam | grep -wf pc2_extreme > pc2_extreme.Q
echo Population Percent Individual > subset_ancestry.txt
cat pc2_extreme.Q | awk '{for (a = 1; a <= 3; a++)print "Pop"a,$a,$4"."$5}' >> subset_ancestry.txt

Rscript -e 'data <- read.table("subset_ancestry.txt",header = TRUE);'\
'library(ggplot2);'\
'bar <- ggplot(data, aes(fill=Population, y=Percent, x=factor(Individual))) +'\
'geom_bar(position="stack", stat="identity") + ggtitle("Admixture ancestries");'\
'bar2 <- bar + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4));'\
'ggsave("subset_ancestries.png", bar2, width = 12, height = 6, dpi = 400);'
```

![](images/subset_ancestries.png)<!-- -->

