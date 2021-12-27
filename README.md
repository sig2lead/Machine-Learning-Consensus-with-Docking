# Machine-Learning-Consensus-with-Docking
1.  Clone Machine-Learning-Consensus-with-Docking repository form Github.
2.  Install requisite R libraries including dplyr, neuralnet, randomForest, caret, and e1071.
3.  Run Sig2Lead connectivity analysis (first tab) and download output “my candidates.csv” file and “lincs candidatess.csv” files and place files in “s2l_files” folder
4.  Run Autodock or other virtual screening program and extract output file and place file in “s2l_files” folder.  The first column should be the compound identifier and the second column should be the median binding energy. 
5.  After placing requisite files in s2l_files folder, run “ml_consensus_with_docking.R” script, after first adjusting "num_compounds" parameter.  The output file will be saved in “results” folder, (you will need to create “results” folder) ordered based on the 80 model Fisher consensus value.
6.  Example files are included in the "s2l_files" folder for EGFR.  These files include "my_candidates.csv", "lincs_candidates.csv", and "docking_file.csv".  The first two files are the output of running the first tab (signature connectivity analysis) of Sig2Lead using EGFR as the gene target.  The last file is an example extracted docking output file for EGFR with compound identifier and median binding energy in the two columns.
