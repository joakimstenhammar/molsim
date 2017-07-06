#!/bin/bash
#*************************************************************
#*      Compare test output with that on Save directory      *
#*************************************************************

sedc='/total cpu time since start/d; /\*\* *version/d; /cpu time/q; /^[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}[ ]\{5\}[0-9]\{2\}:[0-9]\{2\}:[0-9]\{2\}$/d'
i=$1
j=$2
file=`basename $i`
if ! diff -N --ignore-space-change --ignore-blank-lines -q <(sed -e "$sedc" $i) <( sed -e "$sedc" $j) > /dev/null ; then

   echo "differ"
   echo '!!!!!! Start of' $file '!!!!!'
   diff -N --ignore-space-change --ignore-blank-lines <(sed -e "$sedc" $i) <( sed -e "$sedc" $j)
   echo '!!!!!! End of' $file '!!!!!'
   echo  ' '
fi
