#!/usr/bin/env bash

stable="stable.md5sum"
current=".current.md5sum"
dir=$(dirname $1)
file=$(basename $1)
pro=${file%.*}

if [ -f "../version.conf" ]
then
    version="$(cat ../version.conf)"
fi
molsim="$HOME/bin/molsim_ser"
if [ ! -z "$VAR" ]; then
    molsim="$molsim.$version"
fi
if [ "$dir" == "Save" ]; then
    if ! `cmp --silent $stable $current`; then
        echo "The current version is not the stable one. Not creating $1".
        if [ -f $1 ]; then
            echo "Using existing version of $1"; exit 0
        else
            echo "ERROR: $1 must generated with a stable version."; exit 1
        fi
    fi
fi
cd $dir && pwd && $molsim $pro
