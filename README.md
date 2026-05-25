
# Landscape-Aware Hybrid Metaheuristic for Cross-Dock Door Assignment

The software and data in this repository are a snapshot of the software and data used in the research reported in:
*   **Title**: Landscape-Aware Hybrid Metaheuristic for Cross-Dock Door Assignment
*   **Authors**: Wei Yang, Yang Wang, Abraham P. Punnen, and Rafael Martí

---

## Repository Structure

This project is organized into four main directories, containing the complete dataset, source code, execution results, and plotting configurations used in our paper.

### 1. `data/`
This folder contains the complete benchmark datasets for the Cross-dock Door Assignment Problem (CDAP) analyzed in the paper, divided into two distinct subsets:
*   `data/SetA/`: Benchmark instances introduced by Guignard et al. (2012).
*   `data/SetB/`: Benchmark instances introduced by Nassief et al. (2016).

### 2. `results/`
This directory houses the execution logs, statistical output tables, and the data aggregation scripts:
*   `results/SetA/` & `results/SetB/`: Raw text files recording detailed configurations and solution values over **10 independent runs** for each instance (including objectives, solution representation vectors, and runtimes).
*   `results/ablation_study_summary/`: Contains the comparative experimental tables and results for the ablation study.
*   `results/SetA_results_summary.csv` & `results/SetB_results_summary.csv`: Aggregated output sheets containing the Best, Average, and Average Time statistics.
*   `results/StatisticResults.ipynb`: A Python Jupyter Notebook used to parse the raw logs and generate the summary `.csv` sheets in a reproducible order.

### 3. `scripts/`
This directory contains the visualization modules and the actual figures that appear in the paper:
*   **R Language Code**: Scripts configured to read the data tables and plot performance curves, running profiles, and statistical comparisons.
*   **Figures**: Includes all plot output files (such as `Figure3_RunningProfiles.pdf`) featured in the experimental section of the manuscript.

### 4. `src/`
This directory holds the core C++ implementation of the **LA-HM (Landscape-Aware Hybrid Metaheuristic)** algorithm.

---

## Hardware and Software Requirements

### Hardware
*   **CPU**: Multi-core processor (experiments were conducted on Intel(R) Xeon(R) Gold 6226R @ 2.9 GHz).
*   **Memory**: Standard configuration (the algorithm runs on a single thread).

### Software & Library Dependencies
*   **Operating System**: Linux/Unix-based environment is recommended.
*   **C++ Compiler**: GNU GCC 10.2.0 or higher supporting C++11.
*   **Mathematical Programming Solver**: IBM ILOG CPLEX (version V12.9.0 or compatible) with C++ API setup (Concert Technology). CPLEX is required to solve the initial mathematical models $P1$ and $P2$.
*   **Python Environment**: Python 3.x with Jupyter Notebook and `pandas` library to run the statistical aggregation notebook (`StatisticResults.ipynb`) in the `results/` directory.
*   **R Environment**: R environment with common plotting packages (such as `ggplot2`) to execute the visualization scripts in the `scripts/` directory.

---

## Installation and Compilation

To compile the C++ source code, ensure that your CPLEX environment variables (include paths and library paths) are set correctly. Navigate to the `src/` directory and compile the program using GCC:

```bash
g++ -O3 -std=c++11 -DIL_STD \
  -I/opt/ibm/ILOG/CPLEX_Studio129/cplex/include \
  -I/opt/ibm/ILOG/CPLEX_Studio129/concert/include \
  -L/opt/ibm/ILOG/CPLEX_Studio129/cplex/lib/x86-64_linux/static_pic \
  -L/opt/ibm/ILOG/CPLEX_Studio129/concert/lib/x86-64_linux/static_pic \
  -o la_hm main.cpp \
  -lilocplex -lconcert -lcplex -lm -lpthread -ldl
```
*(Note: Please replace `/opt/ibm/ILOG/CPLEX_Studio129` with the actual path of your IBM ILOG CPLEX Optimization Studio installation.)*

---

## Running the Code

The compiled executable accepts the instance name and various algorithmic parameters via command-line flags.

### Command Structure:
```bash
./la_hm -i <instance_name> -a <param_a> -b <param_b> -c <param_d> -d <param_e> -e <param_g>
```

### Argument Descriptions:
*   `-i`: Specifies the target benchmark instance file.
*   `-a`, `-b`, `-c`, `-d`, `-e`: Input parameters for the algorithm (e.g., Maximum iteration in every LA-ITS, SmallDepth, MediumDepth, LargeDepth, pmax).

### Execution Example:
```bash
./la_hm -i SetB_75x20S5.txt -a 500 -b 20000 -c 40000 -d 100000 -e 5
```

---

## Replication and Results Analysis

### 1. Raw Solution Logs
The raw 10-run outputs for each benchmark instance are stored under `results/SetA/` and `results/SetB/`. For details on the inner structure and data format of these solution logs, please refer to the dedicated **`results/README.md`**.

### 2. Result Aggregation
To process the raw output files and compile them into summary CSV tables (`SetA_results_summary.csv` and `SetB_results_summary.csv`) according to the instance sequences defined in `instanceSetA.txt` and `instanceSetB.txt`:

1.  Navigate to the `results/` folder and open the Jupyter Notebook:
    ```bash
    jupyter notebook StatisticResults.ipynb
    ```
2.  Run the notebook cells. The script will calculate the **Best Solution**, **Average Solution**, and **Average Computation Time** for each instance and export the results.

### 3. Visualizations
To reproduce the charts and plots presented in the paper:
1.  Navigate to the `scripts/` directory.
2.  Open and execute the R plotting scripts using your preferred R IDE or terminal command line to re-generate the figures.

---

## License

This project is licensed under the MIT License.