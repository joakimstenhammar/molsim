#!/bin/bash
#*************************************************************
#*      Compare test output with that on Save directory      *
#*************************************************************

diffout='diff.out'
rm -f $diff
sedc='/total cpu time since start/d; /\*\* *version/d; /cpu time/q; /^[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}[ ]\{5\}[0-9]\{2\}:[0-9]\{2\}:[0-9]\{2\}$/d'
i=$1
if ! diff --ignore-space-change --ignore-blank-lines -q <(sed -e "$sedc" out/$i) <( sed -e "$sedc" Save/$i) > /dev/null ; then

   echo "differ"
   echo '!!!!!! Start of' $i '!!!!!'
   diff --ignore-space-change --ignore-blank-lines <(sed -e "$sedc" out/$i) <( sed -e "$sedc" Save/$i)
   echo '!!!!!! End of' $i '!!!!!'
   echo  ' '
fi
