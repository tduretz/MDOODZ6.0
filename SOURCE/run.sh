#!/bin/bash
# print_args.sh
string="$@" 
suffix=".txt"
string=${string%"$suffix"}
echo "Making and running:" $string
make all MODEL=$string OPT=no OMP=no
./Doodzi_$string $string.txt

