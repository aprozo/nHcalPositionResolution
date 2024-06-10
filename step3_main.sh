#!/bin/sh

workdir="/gpfs/mnt/gpfs02/eic/alpro/output/positionResolutionScan/"
mkdir -p $workdir
cd $workdir
angleDir=Phi${GUN_PHI}Theta${GUN_THETA}
mkdir -p $angleDir
cd $angleDir

export GUN_THETA_MIN=$(echo "$GUN_THETA - 0.0001" | bc)
export GUN_THETA_MAX=$(echo "$GUN_THETA + 0.0001" | bc)
export GUN_PHI_MIN=$(echo "$GUN_PHI - 0.0001" | bc)
export GUN_PHI_MAX=$(echo "$GUN_PHI + 0.0001" | bc)
export GUN_MOMENTUM_MIN=$(echo "$GUN_MOMENTUM - 0.00001" | bc)
export GUN_MOMENTUM_MAX=$(echo "$GUN_MOMENTUM + 0.00001" | bc)

source $workdir/epic/install/setup.sh
if [ -f "$DDSIM_FILE" ] && [ "$(stat -c %s "$DDSIM_FILE")" -gt 10000 ]; then
    echo "$DDSIM_FILE exists."
else
    ddsim --compactFile $DETECTOR_PATH/epic_hcal_only.xml --numberOfEvents ${NUMBER_OF_EVENTS} --random.seed $(date +%N) --enableGun --gun.particle ${GUN_PARTICLE} --gun.thetaMin ${GUN_THETA_MIN}*degree --gun.thetaMax ${GUN_THETA_MAX}*degree --gun.phiMin ${GUN_PHI_MIN}*degree --gun.phiMax ${GUN_PHI_MAX}*degree --gun.distribution uniform --gun.momentumMin ${GUN_MOMENTUM_MIN}*GeV --gun.momentumMax ${GUN_MOMENTUM_MAX}*GeV --outputFile ${DDSIM_FILE}
fi
echo "starting eicrecon"

source $workdir/EICrecon/install/bin/eicrecon-this.sh
if [ -f "$EICRECON_FILE" ] && [ "$(stat -c %s "$EICRECON_FILE")" -gt 10000 ]; then
    echo "$EICRECON_FILE exists."
else
    eicrecon $DDSIM_FILE -Ppodio:output_file=${EICRECON_FILE} -Pjana:nevents=${NUMBER_OF_EVENTS} -Ppodio:output_include_collections="MCParticles,HcalEndcapNRawHits,HcalEndcapNRecHits,HcalEndcapNMergedHits,HcalEndcapNClusters,HcalEndcapNTruthClusters,EcalEndcapNRawHits,EcalEndcapNRecHits,EcalEndcapNClusters,EcalEndcapNClusterAssociations,EcalEndcapNTruthClusters,EcalEndcapNTruthClusterAssociations,EcalBarrelScFiClusters,EcalBarrelTruthClusters,EcalBarrelImagingClusters,EcalBarrelClusters,HcalBarrelClusters,HcalBarrelTruthClusters"
    rm -r fieldmaps
fi

echo "starting analysis"
root -l '/gpfs/mnt/gpfs02/eic/alpro/analysis/positionResolutionScan/readHCalRecoReader.C("'${EICRECON_FILE}'" , "outputFill_'${GUN_MOMENTUM}'.root")'
