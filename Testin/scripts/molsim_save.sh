#!/usr/bin/env bash

stable="stable.md5sum"
current=".current.md5sum"
dir=$(dirname $1)
file=$(basename $1)
pro=${file%.*}

molsim="../../Bin/molsim_ser.exe"
if [ "$dir" == "save" ]; then
   if [ -f $dir/$pro.version ]; then
      if `cmp --silent $stable $dir/$pro.version`; then
         echo "This file was created with the currently stable version. Not running molsim again."
         exit 0
      fi
   fi
   if ! `cmp --silent $stable $current`; then
      echo "ERROR: $1 must generated with a stable version."; exit 1
   fi
fi
cd $dir && pwd && $molsim $pro; e=$? && cd -
if [ "$dir" == "save" ]; then
   cat $current > $dir/$pro.version
fi
exit $e
