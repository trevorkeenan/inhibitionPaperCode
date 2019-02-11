# #!/bin/sh

# this is the master script to partition all sites, allowing for inhibition
# this code is structured to run 'batches' 
# to avoid running out of system memory
# each batch uses all cores to run the sites contained in it

# T. Keenan, November 2018 

start=`date +%s`

cd "$HOME/code_batch/code_partitioning/"

batchSize=20		# this is the size of each batch
numSites=213		# total number of sites in FLUXNET

# set the outer loop to run each batch
for yyy in `seq 1 $batchSize $numSites`
do	
	min=$yyy
	max2=$((yyy + batchSize - 1))
	echo '*************************************************NEW BATCH*********************' $min	
	echo '*************************************************NEW BATCH*********************' $min
	echo '*************************************************NEW BATCH*********************' $min	

	# the inner loop runs the sites within each batch
	for xxx in `seq $min 1 $max2`
	do
		indX=$((xxx))
		Rscript ./A1_ScriptPart_main_All_dev2args_wInhib.R $indX & 
	done
 
	wait
 
done

 
end=`date +%s`
runtime=$((end-start))
echo '*******************   Total run-time (s):'
echo "$runtime"
