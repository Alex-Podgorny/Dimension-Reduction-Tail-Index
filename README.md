# Dimension Reduction Tail Index

This supplementary material repository contains code, data, instructions and reproducible results from the article:

> Laurent Gardes and Alex Podgorny. "Dimension reduction for the estimation of the conditional tail-index".  

- Preprint (2024): [hal-04589742](https://hal.science/hal-04589742)

## Project Structure

- **Methods**: Contains R scripts for the estimation of the conditional tail index, including dimensionality reduction methods. This folder includes my proposed method for Conditional Tail Index (CTI) subspace estimation, alongside competitor methods—TDR space (Gardes 2018) and TIREX (Aghbalou et al. 2024).

- **Simulations**: Contains R scripts to generate data for various models, apply the different estimation methods, compute errors, and display results. This folder also includes generated data to ensure reproducibility.

- **Real Data**: Contains the real-world dataset used in the article, along with R scripts for applying my method and displaying results for this dataset.
  
## Observational Data Availability

In the case study application, we use the daily concentration pollutants in Chicagofrom 1987 to 2000. These data were obtained few years ago from the Internet-based Health and Air Pollution Surveillance System (iHAPSS) website but unfortunately they are no longer available. 
____

By Alex Podgorny, IRMA, University of Strasbourg (FR), 2024.   
