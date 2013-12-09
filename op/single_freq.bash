#!/bin/bash

#Initialize Variables
OPDIR="/home/alex/uhd/host/build/examples"
cd ${OPDIR}
BINDIR="/home/alex/sounder"
PROCESS="${BINDIR}/process2.py"
VISUAL="${BINDIR}/visual.py"

rm command_table.txt
rm sounder_comm.fifo

#Create sequence of commands that will be piped to the radar process "sounder"
touch command_table.txt
freq=4499000
txfile="/home/alex/sounder/tx.dat"
rxfile="${OPDIR}/20131126/rx.`TZ=America/Anchorage date +%H%M`.dat"
samprate="250000"
echo "$txfile $rxfile $freq $samprate $samprate" >> command_table.txt
echo "EXIT" >> command_table.txt
cat command_table.txt

#Create pipe, start the USRP radar process, and send commands
mkfifo sounder_comm.fifo
./sounder --args "addr=192.168.10.2" < sounder_comm.fifo &
cat command_table.txt > sounder_comm.fifo

#Wait for the above process to finish
wait $! 

#Remove the command table and the pipe
rm command_table.txt
rm sounder_comm.fifo

#Process the raw data file(s), yielding a numpy array output file
flist=`ls ${OPDIR}/rx.*.dat`
PROCESS="/home/alex/sounder/process3.py"
for f in ${flist}; do
        t=${f%%.dat}
        t=${t##*rx.}
        outfile="rx.${t}"
        echo ${outfile}.dat
        pwd
        if [ -s "${OPDIR}/${outfile}.npy" ]; then continue; fi
        ${PROCESS} -R ${f} -f 400 -n 8192 -O ${outfile}
done

#Remove the raw data file(s)
for f in ${flist}; do
	rm ${f}
done

exit
