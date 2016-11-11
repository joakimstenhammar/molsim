#!/bin/bash - 
#===============================================================================
#
#          FILE: process.diff.out.sh
# 
#         USAGE: ./process.diff.out.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Cornelius Hofzumahaus (CH), hofzumahaus@pc.rwth-aachen.de
#  ORGANIZATION: 
#       CREATED: 11/11/16 10:12:32
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
set +e 
set +x

# ... some declarations
d_out="./out"
d_sav="./Save"
d_bak="./bak"
f_dif="./diff.out"
f_don="./process.done.txt"

# ... Has ./goall.sh been executed?
if [ ! -d $d_out ]; then 
   echo -e "\n You have to run ./goall.sh first! \n\tExiting ..."
   exit
fi

# ... some preliminaries
if [ ! -f $f_don ]; then 
   touch $f_don
   rm -rf $d_bak
else
   echo -e "\n The following files have already been processed:\n"
   echo -e "$( cat $f_don )\n"
fi
mkdir -p $d_bak

# ... Start
echo -e "\nStart of $0\n"

# ... How many files have changed?
nfile=$( grep -e "^differ *$" $f_dif | wc -l )

# ... Process changed files
for ifile in `seq 1 1 $nfile`
do
   # ... file name
   txfile=$( grep --no-group-separator -A 1 -e "^differ *$" $f_dif | head -n $((ifile+ifile)) | tail -n 1 | awk '{ print $4 }')
   if [ $( grep $txfile $f_don | wc -l ) -gt 0 ]; then continue; fi
   lmenu=true
   while $lmenu
   do 
      echo -e "\n! ! ! ! $txfile differs ! ! ! !\n"
      echo -e "How would you like to procede? Choose from\n\ta) vimdiff\n\tb) diff\n\tc) save file\n\td) interrupt\n"
      echo "ENTER CHOICE: ..."
      read -n 1 -r choice
      case $choice in
         "a"|"vimdiff")
            vimdiff $d_out/$txfile $d_sav/$txfile
         ;;
         "b"|"diff") # 'sedc' and diff-command adapted from Pascal
            sedc='/total cpu time since start/d; /\*\* *version/d; /cpu time/q; /^[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}[ ]\{5\}[0-9]\{2\}:[0-9]\{2\}:[0-9]\{2\}$/d'
            diff --ignore-space-change --ignore-blank-lines <(sed -e "$sedc" $d_out/$txfile) <( sed -e "$sedc" $d_sav/$txfile)
         ;;
         "c"|"save")
            cp $d_sav/$txfile $d_bak/${txfile}.bak
            #cp $d_out/$txfile $d_sav/$txfile
            echo -e "\nBackup and overwrite $txfile in $d_sav"
            echo $txfile >> $f_don
            fracdone=$( cat $f_don | wc -l | awk -v nline_awk=$nline '{ printf "%f5.1, $1/nline_awk }')
            echo -e "$fracdone % of changed files have been reviewed.\n"
            lmenu=false
         ;;
         "d"|"interrupt")
            echo "Interrupt script! You can continue from this point at any time."
            exit
         ;;
         *)
            echo "No valid choice. Please enter a valid choice"
         ;;
      esac
   done
done
