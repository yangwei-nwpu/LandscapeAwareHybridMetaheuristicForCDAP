# Detailed Experimental Results of LA-HM on CDAP

This repository contains the detailed experimental results and statistical analysis for the LA-HM algorithm applied to the CDAP (Cross-Docking Assignment Problem). 

---

## Directory Structure

```text
.
├── ablation_study_summary/   # Statistical results of the ablation study experiments
├── SetA/                      # Folder containing detailed results for Set A instances
├── SetB/                      # Folder containing detailed results for Set B instances
├── SetA.tar                   # Compressed archive of the SetA results folder
├── SetB.tar                   # Compressed archive of the SetB results folder
├── instanceSetA.txt           # List of instance names in Set A (defines the summary output order)
├── instanceSetB.txt           # List of instance names in Set B (defines the summary output order)
├── SetA_results_summary.csv   # Statistical summary for Set A instances
├── SetB_results_summary.csv   # Statistical summary for Set B instances
├── detailed_comparison_results.xlsx # Detailed comparison results of LA-HM algorithm and state-of-the-art algorithms on SetA and SetB instances
└── StatisticResults.ipynb     # Jupyter Notebook script for result aggregation
```

---

## Detailed Result Format

The `SetA` and `SetB` folders contain the raw output files for each instance. Every instance is executed **10 times** independently. A typical output block in these files is structured as follows:

```text
O:5174 and computed:5174  and if feasible=1
T:0.000000s
s:1 3 2 1 3 1 0 0 
t:3 1 0 2 0 1 3 2 
***********
```

### Explanation of Fields:
*   **`O`**: The objective value obtained in this run (along with verification of the calculated value and feasibility status `feasible=1`).
*   **`T`**: The computation time taken to find this solution (in seconds).
*   **`s` and `t`**: The representation vectors of the obtained solution for this specific run.
*   **`***********`**: Divider representing the end of one run.

---

## Data Analysis & Summary

The provided script `StatisticResults.ipynb` is used to parse the raw outputs from individual runs and generate a consolidated summary.

### Execution Details:
1. **Input**: The script reads the raw result files from the `SetA/` and `SetB/` directories.
2. **Order of Instances**: The processing order of the instances is strictly determined by the sequence specified in `instanceSetA.txt` and `instanceSetB.txt`.
3. **Metrics Computed**: For each instance (across its 10 runs), the script calculates:
   * **Best Solution**: The best objective value found.
   * **Average Solution**: The mean objective value across the 10 runs.
   * **Average Time**: The mean computation time of the 10 runs.
4. **Output**: The aggregated statistics are exported to `SetA_results_summary.csv` and `SetB_results_summary.csv`.

---

## Ablation Study

The folder `ablation_study_summary` contains additional analytical data and comparative tables representing the performance of different components of the LA-HM algorithm.