# Dimension Reduction Tail Index

This supplementary material repository contains code and reproducible results from the article:

> Laurent Gardes and Alex Podgorny. "Dimension reduction for the estimation of the conditional tail-index".  

- Preprint (2024): [hal-04589742](https://hal.science/hal-04589742)

## Project Structure

This repository is organized into several directories, each corresponding to specific sections and methodologies presented in the accompanying paper.

### **1. Methods**
This folder contains the core scripts implementing the methods used in the study:
- **`local_Hill.R`**: Implements the conditional tail-index estimator (Definition 5).
- **`Minimization_method.R`**: Contains the optimization algorithm described in Section 4.2.
- **`CTI_Estimation.R`**: Implements the CTI (Central Tail-Index) estimator (Definition 6).
- **`Competitors/`**: A subdirectory containing the implementations of competing methods:
  - **`TIREX.R`**: TIREX 1 and 2 proposed in Aghbalou, A., Portier, F., Sabourin, A. and Zhou, C. (2024). Tail Inverse Regression: dimension reduction for prediction of extremes, Bernoulli.
  - **`TDR.R`**: Developed in Gardes, L. (2018). Tail Dimension Reduction for extreme quantile estimation, Extremes.
- **`utils.R`**: Includes utility functions, such as matrix normalization and other auxiliary operations.


### **2. Simulations**
This folder corresponds to the simulation study presented in Section 4 and is organized into three subdirectories:

#### **a. Simulated Data**
- **Purpose**: Covers the data generation process based on the models described in Section 4.3.
- **Contents**:
  - **`Models.R`**: Defines the simulation models used in the study.
  - **`Generating_Data.R`**: Script to generate the simulated datasets.
  - **`Generated Data/`**: Stores all the simulated datasets, ensuring reproducibility.

#### **b. Simulation Results**
- **Purpose**: Includes scripts and results related to the performance analysis in Section 4.4.
- **Contents**:
  - **`Compute.R`**: Performs the simulation computations for the estimation methods.
  - **`Results.R`**: Generates the plots summarizing the results. These are available in the **`Results/`** folder.
  - **Errors file**: A compressed `.zip` file containing all the computed error metrics.

#### **c. Dimension Estimation**
- **Purpose**: Addresses the dimension estimation process described in Section 4.5.
- **Contents**:
  - **`Compute.R`**: Executes the dimension estimation computations.
  - **`Results.R`**: Produces a `Results.RData` file containing the tables presented in the paper.
____

By Alex Podgorny, IRMA, University of Strasbourg (FR), 2024.   
