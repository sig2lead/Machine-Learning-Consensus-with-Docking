# Machine-Learning-Consensus-with-Docking
1.  Clone Machine-Learning-Consensus-with-Docking repository form Github.
2.  Run Sig2Lead connectivity analysis (first tab) and download output “my candidates.csv” file and “lincs candidatess.csv” files and place files in “s2l_files” folder
3.  Run Autodock or other virtual screening program and place output file in “s2l_files” folder.  The first column should be the compound identifier and the second column should be the binding energy. 
4.  After placing requisite files in s2l_files folder, run “ml_consensus_with_docking.R” script, after first adjusting "num_compounds" parameter.  The output file will be saved in “results” folder, (you will need to create “results” folder) ordered based on the 80 model Fisher consensus value.
