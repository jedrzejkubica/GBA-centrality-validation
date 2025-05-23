#!/bin/bash
source ~/pythonvenvs/pyEnv_Multixrank/bin/activate

for PHENO in MMAF NOA BC CC HYPCARD DILCARD OG PCD POF
do
	mkdir $PHENO

	cp default/config.yml $PHENO/.  # make sure the parameters are OK!
	cp -r default/multiplex $PHENO/.  # make sure the interactome TSV file is OK!

	ln -s /home/kubicaj/workspace/GBA-centrality/output/$PHENO/alpha05_d10/ranks_leave_one_out.tsv $PHENO/.
	awk '{print $1}' $PHENO/ranks_leave_one_out.tsv > $PHENO/seeds.txt
	sed -i '1d' $PHENO/seeds.txt

	python run_multixrank.py $PHENO

	python run_leave_one_out.py $PHENO
done
