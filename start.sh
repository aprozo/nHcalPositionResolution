#!/bin/bash

gun_momenta=(0.5 2 10) #gun_momenta
for energy in ${gun_momenta[@]}; do
    ./step0_run.sh $energy
done

~/condor_control.sh # wait for all jobs to finish

# merge different output files (Phi*Theta*) into one
for energy in ${gun_momenta[@]}; do
    workdir=/gpfs/mnt/gpfs02/eic/alpro/output/positionResolutionScan/
    echo "hadd -f -k -j $workdir/output_${energy}.root $workdir/Phi*Theta*/outputFill_${energy}.root" | /eic/u/alpro/eic/eic-shell
done
