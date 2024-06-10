# nHcalPositionResolution

Full script to run is `source run_all.sh`

* `step0_run.sh` - run a single energy (which is an argument) with a chosen $\phi$ and  $\theta$ angles from `tileMap.txt`
it modifies `step1_submit.job` as inserting arguments onto lines 24 and 25

* `step1_submit.job` - condor submit file, one can edit there number of `Queue`'s

* `step2_executable.sh` - is a macro which is executed in a `job` file, it runs `step3_main.sh` in `eic-shell` with some defined variables

* `step3_main.sh` - is the main macro where one :
    *   first, simulation happens `npsim`
    *   `eicrecon`
    *   `root` analysis macro of the eicrecon output