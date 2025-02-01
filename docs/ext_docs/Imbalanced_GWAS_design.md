## Imbalanced GWAS design

A GWAS design with a severely imbalanced cohort (e.g., controls >10x larger than the test cohort) is suboptimal and likely triggers software warnings due to several statistical issues:

***1. Reduced Statistical Power:*** <br>
Power in GWAS is limited by the smaller cohort (test group). Even with a large control group, the ability to detect true associations depends heavily on the test cohort size. The effective sample size is closer to the harmonic mean of the two groups, leading to diminishing returns from an oversized control group.

***2. Biased Variance Estimation:*** <br>
Mixed linear models (MLMs) estimate variance components (e.g., genetic relatedness) using the entire dataset. A dominant control cohort can skew these estimates, leading to inflated/deflated test statistics and spurious associations.

***3. Model Convergence Issues:*** <br>
Extreme imbalance can cause numerical instability during model fitting (e.g., singular matrices, failed convergence). This is especially problematic for MLMs, which rely on iterative algorithms to estimate random effects.

***4. Violation of Homoscedasticity Assumptions:*** <br>
MLMs often assume homogeneous residual variance. If the test cohort has different variance properties (e.g., due to disease status), imbalance exacerbates heteroscedasticity, invalidating results.

***5. Population Stratification Confounding:*** <br>
A disproportionately large control group may introduce hidden population structure if it’s more genetically diverse. While MLMs adjust for relatedness, imbalance can weaken this correction, increasing false positives.
