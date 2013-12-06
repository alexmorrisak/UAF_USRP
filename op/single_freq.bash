#!/bin/bash

OPDIR="/home/alex/uhd/host/build/examples"
BINDIR="/home/alex/sounder"
PROCESS="${BINDIR}/process2.py"
VISUAL="${BINDIR}/visual.py"

cd ${OPDIR}
touch ./freq_table.dat
freq=4096000
txfile="/home/alex/sounder/tx.dat"
starttime=`TZ=America/Anchorage date +%H%M`
rxfile="rx_${starttime}.dat"
samprate="250000"
echo "$txfile $rxfile $samprate $freq" >> freq_table.dat
echo "EXIT" >> ./freq_table.dat
cat freq_table.dat

mkfifo fifo
#./sounder --nsamps 6000000 --args "addr=192.168.10.2" --tx-subdev "A:A" < fifo &
/home/alex/sounder/alexsleep  < fifo &
cat freq_table.dat > fifo 

wait $! 
rm freq_table.dat
rm fifo

cd ${BINDIR}
file=${OPDIR}/${rxfile}
PROCESS="/home/alex/sounder/process3.py"
outfile="rx_${starttime}"
if [ ! -s "/home/alex/sounder/${outfile}.npy" ]; then
  ${PROCESS} -R ${f} -f 400 -n 8192 -O ${outfile}
fi

rm ${rxfile}

exit
