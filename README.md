README – Fiber Photometry and Behavior Analysis (TRPV4 Project)

Author: Tianqi Xu
Date: 2025-07

Project Overview

This repository contains MATLAB codes and experimental datasets used for analyzing fiber photometry (calcium signal) and behavior responses in TRPV4-related ultrasound stimulation experiments. The workflow includes reading raw experiment data, processing and merging behavioral and calcium signals, computing Z-scores, and outputting final analyzed results.


Folder Structure

1. 03_Experiment data

This folder contains the raw experimental data.

Example path:

03_Experiment data\20250113\TX135

This indicates the experiment conducted on January 13, 2025, using the mouse TX135.

2. 01_Matlab Code\03_Github\01_Read data

This folder contains MATLAB scripts for reading and organizing experiment data.

Main script:
FiberBehavirAnalysis_main_v12_250717.m

This is the main file to run the analysis.
Inside the script, fill in the following information:

Experiment date

Mouse ID

Experiment name

Group type

Example:

runFiberBehavior_02_250717('20250113', 'TX135', '120mVpp-2', 'mtTRPV4')

3. Script Configuration
runFiberBehavior_02_250717.m

Line 55: Modify the path of the raw data folder.

Line 65: Modify the path of the output folder where processed results will be saved.

After running this script, the code will automatically read all related data and generate two summary databases:

behaviorData.mat — contains behavioral data.

calciumData.mat — contains calcium fluorescence data.

4. 02_Data analysis Folder

This folder contains scripts for data preprocessing and analysis.

(1) Step10a_CalciumBaseline_01.m & Step10b_MergeMat.m

These scripts merge behaviorData.mat and calciumData.mat into a unified file:

splitData.mat


Each experimental dataset is matched by ID for consistent alignment.

(2) Step11_Zscore.m

Calculates the Z-score of calcium fluorescence signals for normalization.

(3) Step12_Cleanup.m

Extracts baseline data and performs data distribution analysis to compute signal thresholds.

(4) Output

Final behavior and calcium results are summarized and saved as:

outputData.mat


Workflow Summary

1. Read and preprocess raw data → behaviorData.mat, calciumData.mat

2. Merge data → splitData.mat

3. Normalize fluorescence → Z-score

4. Clean up and define thresholds → baseline distribution

5. Output combined results → outputData.mat


Contact

For questions or collaboration, please contact:
Tianqi Xu
Washington University in St. Louis
