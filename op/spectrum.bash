#!/bin/bash

FREQ=${1}
RATE=5000000
FILE="spectrum_${FREQ}_${RATE}.dat"
BINDIR=/home/alex/alex/.pulse/uhd/host/build/examples
BINDIR2=/home/alex/uhd/host/build/examples
$BINDIR/rx_samples_to_file --args addr=192.168.10.2 --wirefmt sc8 --nsamps 262144 --rate ${RATE} --freq ${FREQ} --file ${FILE}
$BINDIR2/spectrum_visual.py --freq ${FREQ} --rate ${RATE} --file ${FILE}
rm ${FILE}
