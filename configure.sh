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
#      REVISION: 2016-10-26 16:39
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

echo -n "checking for FFTW3 using locate .."
if locate -l 1 fftw3.f03 > /dev/null && locate -l 1 libfftw3 > /dev/null ; then # found fftw3
   fftwpaths=`dirname $(locate -b "fftw3.f03") | uniq`
   npath=`echo "$fftwpaths" | wc -l`
   if [ "$npath" -gt "1" ]; then
      echo ""
      echo "Which fftw3 version should be used?"
      select d in $(echo $fftwpaths); do
         if [ -n "$d" ]; then
            echo "$d selected"
            fftwpathfile=$d
            break
         fi
      done
   else
      fftwpathfile=$fftwpaths
   fi
   fftwpath=`dirname $fftwpathfile`
   echo "FFTW_PATH = $fftwpath" >> Src/make.fftwpath
   fftwlibs=`dirname $(locate -r "libfftw3\(.dll\)*\(.a\|.so\)$") | uniq`
   npath=`echo "$fftwlibs" | wc -l`
   if [ "$npath" -gt "1" ]; then
      echo ""
      echo "Which fftw3 lib should be used?"
      select d in $(echo $fftwlibs); do
         if [ -n "$d" ]; then
            echo "$d selected"
            fftwlibfile=$d
            break
         fi
      done
   else
      fftwlibfile=$fftwlibs
   fi
   fftwlib=`dirname $fftwlibfile`
   echo "FFTWLIB = $fftwlib" >> Src/make.fftwpath

elif [ -f "$HOME/.fftw/include/fftw3.f03" ]; then # check for local installation
   echo "FFTW_PATH = $HOME/.fftw" >> Src/make.fftwpath
   fftwlib=`ls -d $HOME/.fftw/lib*`
   echo "FFTWLIB = $fftwlib" >> Src/make.fftwpath
else
   echo "FFTW3 is not found"
   read -e -p "Install under ~/.fftw? (y/n)" -i "n" dofftw
   case ${dofftw:0:1} in
      y|Y )
         FILE="fftw-3.3.4.tar.gz"
         dnfftw=""
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
         if [ -f fftw-*.tar.gz ];
         then
            case ${dnfftw:0:1} in
               y|Y )
                  rm fftw-3.3.4.tar.gz
                  ;;
               * )
                  echo "not removing tar file"
            esac
         fi
         cd $HOME/.fftw/fftw*
         pwd
         ./configure --prefix="$HOME/.fftw"
         make
         make install
         cd $curdir
         echo ""
         echo "FFTW installed"
         echo ""
         echo "FFTW_PATH = $HOME/.fftw" >> Src/make.fftwpath
         fftwlib=`ls -d $HOME/.fftw/lib*`
         echo "FFTWLIB = $fftwlib" >> Src/make.fftwpath
         ;;
      * )
         echo "Please provide the path to the directory of the fftw3.f03 file:"
         read -i "/usr/local/" fftwpath1
         fftwpath=`dirname $fftwpath1`
         echo "FFTW_PATH = $fftwpath" >> Src/make.fftwpath

         echo "Please provide the path to the directory of the fftw3 libary (libfftw3):"
         read -i "$fftwpath/lib" fftwlib1
         fftwlib=`dirname $fftwpath1`
         echo "FFTWLIB = $fftwlib" >> Src/make.fftwpath
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


#creating version.conf where user parameters are stored
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
