## GCTA-GREML Analysis

The narrow-sense heritability (h<sup>2</sup>) is defined as the proportion of phenotypic variance due to additive genetic variance. GREML (Genomic-relatedness-based restricted maximum likelihood) is a method used to estimate the proportion of variance in a phenotype that can be explained by SNPs in a population. In other words, **GREML calculates the SNP-based heritability**.

GREML utilizes the genetic relationships between individuals in a sample population, captured in a genomic relationship matrix (GRM), to *partition the phenotypic variance* into its genetic and residual components. By fitting a linear mixed model, GREML estimates the proportion of variance attributable to all SNPs.

Here is a deeper dive into linear mixed models (LMM) and how they are used in the context of this genomic analysis:

***1. Basics of linear models:*** <br>
In simple linear regression, the relationship between a dependent variable (phenotype) and one or more independent variables (often genotypic data) is modeled. The model is described as:

Y=Xβ+ϵ

Where:

*  Y is the vector of phenotypic values.
*  X is the design matrix for the fixed effects.
*  β is a vector of fixed effect coefficients.
*  ϵ is the vector of random errors.

<br>

***2. Mixed models add a random component:*** <br>
In mixed models, besides fixed effects (like the β in the linear model), there are random effects which capture variability in the data due to certain groups or clusters.

Y=Xβ+Zu+ϵ

Where:

*  Z is the design matrix for the random effects.
*  u is the vector of random effects.

<br>

***3. Application to genomic data:*** <br>
In the context of GREML and genomic studies:

*  The fixed effects (β) might capture known factors influencing the phenotype, like age or gender.
*  The random effects (u) capture the genetic effects of all SNPs. These are assumed to follow a multivariate normal distribution with mean zero and variance-covariance matrix Gσ<sup>2</sup><sub>g</sub>, where G is the genomic relationship matrix (GRM) and σ<sup>2</sup><sub>g</sub> is the genetic variance.
*  The residuals ϵ are also assumed to be normally distributed with mean zero and a variance-covariance matrix Iσ<sup>2</sup><sub>e</sub>, where I is the identity matrix and σ<sup>2</sup><sub>e</sub> is the error variance.

<br>

***4. Estimating the variance components:*** <br>
The key goal in GREML is to estimate the σ<sup>2</sup><sub>g</sub> (genetic variance) and σ<sup>2</sup><sub>e</sub> (residual variance). From these, the SNP-based heritability can be computed:

h<sup>2</sup> = σ<sup>2</sup><sub>g</sub> / (σ<sup>2</sup><sub>g</sub> + σ<sup>2</sup><sub>e</sub>)

This quantifies the proportion of the phenotypic variance explained by the SNPs.

<br>

***5. Restricted Maximum Likelihood (REML):*** <br>
The REML approach is used to estimate the variance components in the model. REML provides unbiased estimates of variance components under a wide range of situations.

In summary, by using a linear mixed model, the GREML approach leverages both fixed (known factors) and random (genetic effects of SNPs) components to partition and quantify the variability in the observed phenotype. This enables researchers to pinpoint how much of the phenotypic variance can be attributed to genetics as captured by SNPs.

<br>


***6. Statistical verification:*** <br>
For GCTA to assess the importance of the SNPs' effects to the model, it fits and compare 2 models:

a) ***Full Model*** that includes both the fixed effects and the random effects. <br>
b) ***Reduced Model*** that only includes the fixed effects and omits the random effects of SNPs. <br>

GCTA calculates the ***Log Likelihoods*** for both models. A higher log likelihood indicates a better fit of the model to the data. The difference in log likelihoods between the full and reduced models provides a test statistic for the significance of the random effects (SNP effects). Essentially, if the full model (which includes SNP effects) has a significantly better fit to the data than the reduced model (which excludes SNP effects), it implies that the SNPs are explaining a significant portion of the variance in the phenotype.

Formally, the difference in log likelihoods between the two models, multiplied by two, approximately follows a chi-squared distribution. This can be used for hypothesis testing:
2 (LogLikelihood<sub>full</sub> − LogLikelihood<sub>reduced</sub>).




More details about linear mixed models can be found [here](https://www.youtube.com/playlist?list=PL8F480DgtpW9_IT7xN1XeRF_dglZmK0nM).
