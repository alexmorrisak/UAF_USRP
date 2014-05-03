#!/bin/bash 

#Initialize Variables
TOP=/home/alex/UAF_USRP
OPDIR=${TOP}/op
cd ${OPDIR}
BINDIR=${TOP}/process
DATADIR=/home/alex/ionosonde_data
SOUNDER="${TOP}/uhd/sounder/sounder"
PROCESS="${BINDIR}/process.py"
VISUAL="${BINDIR}/visual.py"
DATE=`date +%Y%m%d`
TIME=`date +%H%M`

mkdir -p ${DATADIR}/${DATE}/${TIME}

rm command_table.txt
rm sounder_comm.fifo
rm process_in.fifo

touch command_table.txt
for i in {1..50}; do
        freq=`expr 2400000 \+ $i \* 100000`
	freq=`printf %.8d ${freq}`
        txfile="/home/alex/UAF_USRP/process/tx.dat"
	fsize=`stat -c%s $txfile`
	nsamples=`expr $fsize \/ 4`
	rxfile="${DATADIR}/${DATE}/${TIME}/rx.${TIME}.${freq}.dat"
        samprate="250000"
        echo "$txfile $rxfile $nsamples $freq $samprate $samprate" >> command_table.txt
done
echo "EXIT" >> command_table.txt
cat ./command_table.txt

mkfifo process_in.fifo
${PROCESS} -n 1024 -p MLS_31 < process_in.fifo &

#Create pipe, start the USRP radar process, and send commands
mkfifo sounder_comm.fifo
echo "made fifo"
${SOUNDER} > process_in.fifo < sounder_comm.fifo &
echo "started sounder"
cat command_table.txt > sounder_comm.fifo

#Wait for the above process to finish
wait %1
wait %2

#Remove the command table and the pipe
rm command_table.txt
rm sounder_comm.fifo
rm process_in.fifo

exit
