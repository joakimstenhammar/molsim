#!/bin/bash
#*************************************************************
#*      Compare test output with that on Save directory      *
#*************************************************************
 
diffout='diff.out'
rm -f $diffout
for i in `cat todiff.txt`
do
#   echo '!!!!!! Start of' $i '!!!!!'
   ./diff1.sh $i >> $diffout
done
 
less $diffout
