#!/bin/bash
GUN_MOMENTUM=${1}
workdir="/gpfs/mnt/gpfs02/eic/alpro/output/positionResolutionScan"
current_dir=$(pwd)

mkdir -p $workdir
cd $workdir

################################################################################################
#  make compiled EPIC and EICrecon
################################################################################################

# if [ -d "epic" ]; then
#   echo "EPIC already exists."
# else
#   cp -r /gpfs/mnt/gpfs02/eic/alpro/output/epic .
#   cd epic
#   sed -i '59s/grid_size_x="[^"]*"/grid_size_x="'$GRID_SIZE'0 * mm"/' compact/hcal/backward.xml
#   sed -i '59s/grid_size_y ="[^"]*"/grid_size_y="'$GRID_SIZE'0 * mm"/' compact/hcal/backward.xml

#   cat <<EOF | /eic/u/alpro/eic/eic-shell
#   echo "Compiling EPIC"
#   cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install
#   cmake --build build -j8 -- install
#   source install/setup.sh
#   cd ..
#   git clone https://github.com/eic/EICrecon
#   cd EICrecon
#   cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install
#   cmake --build build
#   cmake --install build
# EOF
# fi

################################################################################################

cd $current_dir
# Define the input file
# input_file="tileMap.txt"
input_file="tileMapTest.txt"
# Define the base job configuration file
base_config="step1_submit.job"

echo "reading from $input_file"

# Check if the input file exists
if [ ! -e "$input_file" ]; then
  echo "Input file '$input_file' does not exist."
  exit 1
fi
mkdir -p $workdir/log

################################################################################################
# Loop through each line in the input file to run over many different PHI and THETA values
################################################################################################
while read -r line; do
  # Split the line into two variables, assuming they are separated by a space
  PHI=$(echo "$line" | cut -d' ' -f1)
  THETA=$(echo "$line" | cut -d' ' -f2)

  echo "PHI: $PHI"
  echo "THETA: $THETA"

  # Define an array of argument sets
  arg_set=("${PHI} ${THETA} ${GUN_MOMENTUM}")
  # Create a new job configuration file by appending arguments to the base file
  config_file="config_$(echo $arg_set | tr ' ' '_').job"
  cp "$base_config" "$config_file"

  line="Outputdir               = ${workdir}"
  sed -i "24 a ${line}" $config_file

  line=" arguments =  ${arg_set} \n"
  sed -i "25 a ${line}" $config_file

  # Submit the job
  echo "Submitting job with arguments: $arg_set"
  condor_submit "$config_file"
  # remove the temporary config file
  rm "$config_file"

done <"$input_file"
