# - Variables of interest:
#   - Number of CpGs in region $K$
#   - Number of cells $N$
#   - Level of CpG density
#     - Rich: CpG-CpG distance ~ $\text{Poisson}(\lambda = 10)$
#     - Moderate:  CpG-CpG distance ~ $\text{Poisson}(\lambda = 250)$
#     - Poor: CpG-CpG distance ~ $\text{Poisson}(\lambda = 1,000)$
#   - Methylation info missing rate $\gamma$ (hence average across-cell coverage should be $N\gamma$ )
#   - Noise level $\sigma$
#   - Prevalence $\pi$ (only for two grouping case)