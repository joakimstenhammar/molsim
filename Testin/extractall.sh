#!/bin/bash - 
#===============================================================================
#
#          FILE: extractall.sh
# 
#         USAGE: ./extractall.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Pascal Hebbeker (PH), pascal.hebbeker@gmail.com
#  ORGANIZATION: 
#       CREATED: 2016-08-15 10:12
#      REVISION: 2016-08-15 14:03
#===============================================================================

set -o nounset                              # Treat unset variables as an error

expand () {
   xz -kvd $1.tar.xz
   tar -xvf $1.tar
   rm $1.tar
}

echo "expanding in"
expand in
echo "expanding Save"
expand Save
