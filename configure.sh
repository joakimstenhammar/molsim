#!/bin/bash - 
#===============================================================================
#
#          FILE: configure.sh
# 
#         USAGE: ./configure.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (), 
#  ORGANIZATION: 
#       CREATED: 28/07/2016 09:09
#      REVISION: 2016-10-26 12:53
#===============================================================================

set -o nounset                              # Treat unset variables as an error
set -e

i=ifort
echo -n "Checking ifort ..."
command -v ifort >/dev/null 2>&1 || {
   echo "no";
   echo -n "Checking gfortran ..."
   command -v gfortran >/dev/null 2>&1 || {
      echo "Error: fortran compiler required for MOLSIM"
      exit 1;
   }
   echo "yes"
   read -e -p "Use gfortran instead? " -i "y" dogfortran
   case ${dogfortran:0:1} in
          y|Y )
            echo -n "Enabling gfortran in the makefile ..."
            echo "ARCH = LOCAL_GFORTRAN" >> Src/make.arch
         ;;
         * )
            echo "Error: fortran compiler required for MOLSIM"
            exit 1;
      esac
}
echo "yes"

#progs=("")
#for i in "${progs[@]}"; do
   #echo -n "Checking $i ..."
   #command -v $i >/dev/null 2>&1 || { echo >&2 "I require $i but it's not installed.  Aborting."; exit 1; }
   #echo "yes"
#done

echo -n "checking FFTW3 .."
if [ ! -f "$HOME/.fftw/include/fftw3.f03" ]; then
   echo "FFTW3 is not installed under ~/.fftw"
   read -e -p "Install? " -i "n" dofftw
   case ${dofftw:0:1} in
       y|Y )

      FILE="fftw-3.3.4.tar.gz"

      if [ ! -f $FILE ];
      then
         read -e -p "$FILE not found. Download from fftw.org? (requires wget) " -i "n" dnfftw
         case ${dnfftw:0:1} in
             y|Y )
             wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz
             ;;
             * )
            read -e -p "Provide path to file: " FILE
         esac
      fi

      if [ ! -f $FILE ];
      then
         echo "ERROR: $FILE not found."
         exit 1
      fi

      echo "Using $FILE"
      echo "Installing in $HOME/.fftw"

      mkdir -p $HOME/.fftw
      curdir=$PWD
      tar xfv $FILE -C $HOME/.fftw
      cd $HOME/.fftw/fftw*
      pwd
      ./configure --prefix="$HOME/.fftw"
      make
      make install
      cd $curdir
      ;;
      * )
         echo "Error: FFTW Required for MOLSIM"
   esac
fi
echo "yes"

echo -n "Checking ~/bin ..."
if [[ ! ":$PATH:" == *":$HOME/bin:"* ]]; then
   echo "Setting up ~/bin"
   mkdir -p $HOME/bin
         
   #Add ~/bin to PATH
   export PATH=".:$HOME/bin:$PATH"
   rc=${SHELL#*/bin/}rc
   echo 'PATH=".:$HOME/bin:$PATH"' >> "$HOME/.$rc"
   echo 'export PATH' >> "$HOME/.$rc"
fi
echo "yes"


#creating pdot.conf where user parameters are stored
conffile="version.conf"

if [[ -e "$conffile" ]]; then
   read -e -p "$conffile exists. Overwrite? " -i "n" doconf
else
   doconf="y"
fi

case ${doconf:0:1} in
    y|Y )
      read -e -p "Name of the version? " -i "" ver
      echo $ver > $conffile
    ;;
    * )
        echo "not changing $conffile"
    ;;
esac
