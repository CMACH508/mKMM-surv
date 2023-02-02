# Kernel Mean Matching with Mahalanobis Distance for Causal Inference of Time-to-event Outcome

A number of existing causal effect estimation approaches for time-to-event outcomes in observational studies employ the propensity score-based strategies to balance covariates in treatment groups and reduce confounding effects. The difficulty of causal inference using propensity score lies in the issue that correctly specified propensity score model is required to obtain causal effects, otherwise, substantial bias may be introduced. In this paper, we develop a nonparametric weighting approach, which balances the covariates between treatment groups by Kernel Mean Matching (KMM) with Mahalanobis distance kernels. The proposed estimator adapts KMM procedure to causal effect estimation for time-to-event outcomes, rather than its original applications on empirical risk for classification or regression. Given universal kernels, KMM procedure leads to the balance of covariate distributions. For time-to-event outcome which has highly nonlinear relations with covariates, balancing covariates at distribution level is necessary rather than finite order moments. Furthermore, in consideration of practical applications, unlike the Euclidean distance, the Mahalanobis distance is essentially suitable for the anisotropic real-world data with extremely different scales and highly correlated dimensions, so that we can make the best use of the optimization power. We theoretically proved the consistency of the proposed estimator. Simulations and real-world applications illustrate the effectiveness of the proposed method.

## Requirements

- R 4.0.5
- survival
- kernlab
- quadprog
- Matrix
- MASS
- CBPS
- parallel
- sjstats
- car
- survtmle

## Run the experiments

```shell
Rscript survival_experiment.R
```

## Citation


Please cite our work if you find our code/paper is useful to your work.


```
@INPROCEEDINGS{9995021,

  author={Ma, Qin and Zeng, Lin and Tu, Shikui and Xu, Lei},

  booktitle={2022 IEEE International Conference on Bioinformatics and Biomedicine (BIBM)}, 

  title={Kernel Mean Matching with Mahalanobis Distance for Causal Inference of Time-to-event Outcome}, 

  year={2022},

  volume={},

  number={},

  pages={509-514},

  doi={10.1109/BIBM55620.2022.9995021}}

```


