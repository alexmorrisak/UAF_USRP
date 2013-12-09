#!/bin/bash

# Not a working program.  Demonstrates perhaps how a swept-frequency
# ionosonde might operate.

touch ./freq_table.dat
for i in {1..10}; do
	freq=`expr $i \* 1000000`
	txfile="/home/alex/gnuradio/gnuradio/gr-radar/tx.dat"
	rxfile="rx_`printf "%06d" $freq`.dat"
	samprate="500000"
	echo "$txfile $rxfile $samprate $freq" >> freq_table.dat
done
echo "EXIT" >> ./freq_table.dat
#cat freq_table.dat
#mkfifo fifo
#./sounder --tx-rate 500000 --rx-rate 500000 --freq 1000000 --nsamps 1500000 < fifo &
#cat freq_table.dat > fifo &
rm freq_table.dat
rm fifo
exit
