# Case 3: Team 06

This repository contains the source code, data, and report for our project. Our project focuses on analyzing and visualizing the data using various statistical methods and techniques.

## Repository Structure

Here is a brief overview of the files in this repository:

- `all_code.R`: This file contains all the R code used for data preprocessing, analysis, and visualization. You can find the source code for all the functions, models, and plots used in our project in this file.

- `report_submission1.pdf`: This is the PDF version of our report, which includes a comprehensive explanation of our methodology, results, and conclusions. The report provides detailed insights into our data analysis and findings, as well as recommendations for further research.

- `report_submission1.qmd`: This file is the R Markdown source code for our report. The Rmd file contains all the text, code, and output used to generate the PDF version of the report (`report_submission1.pdf`). You can use this file to reproduce our report or modify it to suit your needs.

## Getting Started

To get started with our project, you can follow these steps:

1. Clone this repository to your local machine.
2. Open the `all_code.R` file in your preferred R environment (e.g., RStudio) and run the code to perform the data analysis and generate the plots.
3. If you wish to regenerate or modify the report, open the `report_submission1.Rmd` file in your R environment and use the "Knit" button to compile the

## Introduction to Project

Glaucoma, a group of eye diseases damaging the optic nerve, can lead to vision loss or irreversible blindness
if left untreated (Quigley & Broman, 2006). It affects nearly 80 million people as of 2023 (Tham et al., 2014).
The disease often progresses slowly, without noticeable symptoms, until irreversible vision loss occurs, making
early detection and monitoring crucial for effective treatment (Kass et al., 2002).
Patients routinely undergo visual field examinations to monitor glaucoma status, providing a functional assessment of vision across the field (Bengtsson & Heijl, 2005). These longitudinal series of visual fields are used
by clinicians to determine disease progression rates. However, the data contains complex spatial and temporal
dependencies that must be considered when determining progression.

In this case study, we investigate early onset glaucoma patients in the dataset and predict glaucoma progression.
We explore whether a spatial model, specifically an areal model, improves prediction over a non-spatial model,
such as simple linear regression. Specifically, we look to train the models on the first 5 visual fields and see if
the spatial or non-spatial model better predicts the future visual fields. This research question holds significant
clinical importance as it can potentially enhance disease progression predictions, leading to improved treatment
options and patient outcomes (Yousefi et al., 2015).

We hypothesize that, among early onset glaucoma patients, spatial modeling using areal models will yield more
accurate predictions of disease progression than non-spatial modeling approaches like simple linear regression.
This is based on the spatial nature of visual field data in glaucoma patients, as the disease affects different
visual field regions in distinct patterns. Areal models, a type of spatial model, account for spatial dependencies
and topographical relationships between visual field regions, revealing insights about disease progression. On
the other hand, non-spatial models do not consider these spatial relationships and may fail to capture essential
patterns.

In general, our results could improve clinical decision-making and care quality for this patient population

To read our results, view Final_Report.pdf. If you would like to look at the code interwoven with the written results,
view Final_Report_With_Code.qmd.
