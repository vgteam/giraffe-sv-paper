# real_read_experiments

These scripts are used for mapping real read HG002 and HG003 GIAB samples against various linear and graph-based alignment algorithms.

To run them, you will need a machine with the `aws` command configured.

You will also need `Docker` installed and runable within the same machine.

For the DRAGEN script, you will need a machine with the Illumina DRAGEN module installed and runable.

Files are output locally on the machine.
`bwamem_map.sh` outputs files to the `${HOME}/run_bwamem_mapping` directory.
`dragen_map.sh` outputs files to the `${HOME}/run_dragen_mapping` directory.
`vg_map.sh` outputs files to the `${HOME}/run_vg_map_mapping` directory.
`giraffe_map.sh` outputs files to the `${HOME}/run_giraffe_mapping` directory.

## Running the scripts

To run BWA-MEM on HG002 and HG003 150bp and 250bp reads:
`./bwamem_map.sh`

To run DRAGEN on HG002 and HG003 150bp and 250bp reads:
`./dragen_map.sh`

To run VG MAP on HG002 and HG003 150bp and 250bp reads:
`./vg_map.sh`

To run VG Giraffe, VG Giraffe in FAST mode, and VG Giraffe on Primary graphs for HG002 and HG003 150bp and 250bp reads:
`./giraffe_map.sh`


