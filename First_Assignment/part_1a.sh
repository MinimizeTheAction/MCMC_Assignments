#!/bin/bash 

make studentt_sampler

# arguments to pass to ./studentt_sampler
no_samples=5000
nu=3.0
mu=1.0
sigma=1.0
data_files=("data_file_0p01.txt" "data_file_0p1.txt" "data_file_1.txt" "data_file_10.txt")
proposal_sigma=(0.01 0.1 1.0 10.0)

# arguments to pass to post_process.py
figure_files=("chain_0p01.pdf" "chain_0p1.pdf" "chain_1.pdf" "chain_10.pdf")


for i in `seq 0 3`;
do
	./studentt_sampler $no_samples $nu $mu $sigma "${data_files[$i]}" "${proposal_sigma[$i]}"
	./post_process.py "${data_files[$i]}" "${figure_files[$i]}" &

done

wait
