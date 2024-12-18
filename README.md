# Systematic Comparison of Mendelian Randomization Methods

## Project Overview
This repository hosts code for a project designed to evaluate and compare a unique set of Mendelian Randomization (MR) methods and Polygenic Risk Score (PRS) methods across diverse scenarios, including both simulations and real-world data. The project examines common and challenging conditions such as weak instruments, linkage disequilibrium, and pleiotropy.

## Noted：Sample Code is a Subset of Whole Project
In practical project implementation, additional code and more complex processes are involved. Displaying everything would make this GitHub project overly complex and less intuitive. Therefore, this repository includes only a selected subset of scripts intended to illustrate the overall workflow.The full codebase would be considerably larger than the example provided here.

### Goals

- Simulate data.
- Calculate the estimated MR effects.
- Compare the performance of various MR and PRS methods under different conditions.

### Code Flowchart
![Code Flowchart](./Figs/Code_flowchart.jpg)


### Methods to be Assessed

The MR methods evaluated include:
- `mr_sisVIVE`
- `mr_2SLS`
- `mr_ivw`
- `mr_median`
- `mr_mode`
- `mr_mix`
- `mr_egger`
- `mr_raps`

The PRS Methods evaluated include:
- `PRS_PT`
- `PRS_Oracle`
- `PRS_LDPred`

The calculated PRS is used as the independent instrumental variable (IV) for subsequent Mendelian Randomization analyses.

## Software and Key Packages

**R Version**: R/4.3.0-foss-2020b  
**Python Version**: Python/3.8.6-GCCcore-10.2.0  

### R Packages
- `CorBin[1]`
- `MASS`
- `MendelianRandomization`
- `sisVIVE`
- `RobustIV`
- `ivreg`
- `genio`


### Python Packages
- `LDpred[2]`
- `SciPy-bundle/2020.11-foss-2020b`


## Package References
[1] Jiang, W., Song, S., Hou, L., Zhao, H., 2020. A Set of Efficient Methods to Generate High-Dimensional Binary Data With Specified Correlation Structures. The American Statistician, 75(3), pp.310–322. [(R CorBin package)](https://cran.r-project.org/web/packages/CorBin/index.html)


[2] Vilhjálmsson, B.J., Yang, J., Finucane, H.K., et al, 2015. Modeling Linkage Disequilibrium Increases Accuracy of Polygenic Risk Scores. The American Journal of Human Genetics, 97(4), pp.576-592..
[(LDpred GitHub)](https://github.com/bvilhjal/ldpred)



