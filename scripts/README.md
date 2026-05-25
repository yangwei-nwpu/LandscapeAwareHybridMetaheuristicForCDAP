# Reproducible Figures for LA-HM Evaluation



This repository contains the R code and datasets required to reproduce the key figures presented in the paper and its online supplement. 



All visualization procedures are implemented in the file `Figure.R`.



---



## Prerequisites



To execute the scripts, ensure you have R installed (the code was developed and tested using \*\*R version 3.6.1\*\*). The following packages are required:



* `ggplot2` (v3.2.1 or later)

* `reshape2` (v1.4.3 or later)

* `ggsci` (v2.9 or later)

* `dplyr` (v0.8.3 or later)

* `cowplot` (v1.0.0 or later)

* `grid` (built-in R library)



You can install the necessary packages using the R command below:



```R

install.packages(c("ggplot2", "reshape2", "ggsci", "dplyr", "cowplot"))

```



---



## Repository Structure



To run the scripts without path modifications, please maintain the directory structure as shown below:



```text

.

├── Figure.R                        # R script containing visualization codes

├── SetBgapLNSVerifyRev.csv         # Input data for the main boxplot (Figure OS1)

├── combined\_csv/                   # Directory containing algorithm performance CSV files

│   ├── instance\_1.csv              # Example data file for the running profiles

│   ├── instance\_2.csv

│   └── ...

└── README.md                       # Documentation file

```



---



## Details of Figures



### 1. Running Profiles of LA-HM with its Variants

* **Paper Reference**: Section 4.3, Figure 3 ("Running profiles of LA-HM with its four variants")

* **Input Data**: Located in the `combined\_csv/` directory. Each file contains convergence data over time for different algorithm variants.

* **Output**: Individual line charts saved as PDF files (with the suffix `\_JCO.pdf`) inside the `combined\_csv/` folder.



### 2. Boxplot of Percent Gaps

* **Paper Reference**: Online Supplement Section 2, Figure OS1 ("Boxplot of percent gaps among LA-HM, LA-ITS and MLNS for best and average objective values")

* **Input Data**: `SetBgapLNSVerifyRev.csv` (located in the root directory).

* **Output**: `SetBGapLNSVerigy.pdf` (saved in the root directory). This figure contains a main boxplot with an embedded inset plot to focus on specific comparisons.



---



## Execution Instructions



1\. Set the working directory in your R session to the repository's root folder:

  ```R

  setwd("path/to/your/repository")

  ```

2\. Run the code from the command line or source the script directly in your R console:

  ```R

  source("Figure.R")

  ```

3\. Alternatively, you can open Figure.R and execute the code blocks sequentially for Part 1 (Running Profiles) and Part 2 (Boxplot with Inset).

