# MLTSKC-FS

Repository associated with the manuscript:

**“Enhancing Multi-label Classification by a Choquet Integral Based Generalization of the Multi-label Takagi–Sugeno–Kang Fuzzy System Model”**

submitted to **Applied Soft Computing**.

## Authors

- Karina Condori
- Julian Suarez
- Giancarlo Lucca
- Qiongdan Lou
- Zhaohong Deng
- Humberto Bustince
- Graçaliz P. Dimuro

## About the project

This repository contains the data, source code, and supplementary material used in the study of the **Multi-Label Takagi–Sugeno–Kang Choquet Fuzzy System (ML-TSKC FS)**.

The proposed model generalizes the antecedent rule-activation mechanism of the original ML-TSK FS by replacing the product-based aggregation of antecedent membership degrees with a discrete Choquet-integral-based aggregation mechanism. Five fuzzy-measure configurations are investigated:

- Uniform
- Relative
- Product
- Power
- Weighted

The repository is intended to support the reproducibility of the experiments reported in the manuscript.

## Repository structure

- **Data/**  
  Contains the datasets used in the experimental study.

- **Source code/**  
  Contains the MATLAB (`.m`) source code used to implement and evaluate ML-TSKC FS, including the Choquet-integral-based aggregation procedure.

- **Additional Tables/**  
  Contains supplementary tables and additional experimental material related to the manuscript.

## Experimental protocol

The experiments use 12 multi-label benchmark datasets and four evaluation metrics:

- Average Precision (AP)
- Hamming Loss (HL)
- Ranking Loss (RL)
- Coverage (CV)

The adjustable hyperparameters are selected by grid search combined with five-fold cross-validation, following the experimental protocol adopted in the original ML-TSK FS study. The fuzzy measures themselves are not learned during consequent optimization.

## Relation to ML-TSK FS

ML-TSKC FS is developed as a generalization of the Multi-Label Takagi–Sugeno–Kang Fuzzy System proposed by Lou et al.:

> Q. Lou, Z. Deng, Z. Xiao, K.-S. Choi, and S. Wang,  
> “Multilabel Takagi-Sugeno-Kang fuzzy system,”  
> *IEEE Transactions on Fuzzy Systems*, vol. 30, pp. 3410–3425, 2022.

The present work retains the first-order multi-label TSK structure while modifying the antecedent aggregation mechanism through the Choquet integral.

## Reproducibility

The repository provides the material required to reproduce the experiments described in the manuscript. Please refer to the source-code files and supplementary material for the corresponding implementation and experimental settings.

## Citation

If you use this repository or the ML-TSKC FS implementation in your research, please cite the associated manuscript. Full bibliographic information will be added after publication.

## License

Please refer to the repository license file, if provided, for terms of use and redistribution.
