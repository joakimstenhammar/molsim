#!/usr/bin/env bash
set -o nounset                              # Treat unset variables as an error
set -e

echo -n "Checking ifort ..."
if command -v ifort >/dev/null 2>&1; then
   echo "yes"
   lifort=true
else
   echo "no"
   lifort=false
fi
echo -n "Checking gfortran ..."
if command -v gfortran >/dev/null 2>&1; then
   echo "yes"
   lgfortran=true
else
   echo "no"
   lgfortran=false
fi
echo -n "Checking mpiifort ..."
if command -v mpiifort >/dev/null 2>&1; then
   echo "yes"
   lintelmpi=true
   intelmpiFC=$( mpiifort -show | { read -a array ; echo ${array[0]} ; })
   echo "mpiifort uses $intelmpiFC"
else
   echo "no"
   lintelmpi=false
   intelmpiFC=""
fi
echo -n "Checking mpifort ..."
if command -v mpifort >/dev/null 2>&1; then
   echo "yes"
   lopenmpi=true
   openmpiFC=$(mpifort -show | { read -a array ; echo ${array[0]} ; })
   echo "mpifort uses $openmpiFC"
else
   echo "no"
   lopenmpi=false
   openmpiFC=""
fi

FC="none"
MPIFC="none"
if [ $lifort = true ] && [ "$intelmpiFC" = "ifort" ]; then
   echo "Found intelmpi and ifort"
   FC="ifort"
   MPIFC="mpiifort"
elif [ $lifort = true ] && [ "$openmpiFC" = "ifort" ]; then
   echo "Found openmpi and ifort"
   FC="ifort"
   MPIFC="mpifort"
elif [ $lgfortran = true ] && [ "$openmpiFC" = "gfortran" ]; then
   echo "Found openmpi and gfortran"
   FC="gfortran"
   MPIFC="mpifort"
elif [ $lgfortran = true ] && [ "$intelmpiFC" = "gfortran" ]; then
   echo "Found intelmpi and gfortran"
   FC="gfortran"
   MPIFC="mpiifort"
else
   echo "Warning: Could not detect a working mpi compiler combination. Unless the mpi compiler is adapted, the compilation of the parallel version will fail"
   if [ $lifort = true ]; then
      echo "Using ifort"
      FC="ifort"
   elif [ $lgfortran = true ]; then
      echo "Using gfortran"
      FC="gfortran"
   else
      echo "Warning: Automatic detection of the fortran compiler failed"
   fi
fi

setcomp=false
while [ $setcomp = false ]; do
   read -e -p "Use $FC as a compiler and $MPIFC as mpi compiler? " -i "y" docomp
   case ${docomp:0:1} in
      y|Y )
         setcomp=true
         ;;
      * )
         read -e -p "Which compiler to use? " -i "$FC" FC
         read -e -p "Which mpi compiler to use? " -i "$MPIFC" MPIFC
   esac
done
if [ "$FC" = "gfortran" ]; then
   echo "ARCH = LOCAL_GFORTRAN" > Src/make.arch
elif [ "$FC" = "ifort" ]; then
   echo "ARCH = LOCAL_INTEL" > Src/make.arch
else
   echo "Compiler $FC is not supported by molsim"
   exit 1
fi
echo "MPIFC = $MPIFC" >> Src/make.arch

fftwpath="other"
fftwlib="other"
echo -n "checking for FFTW3 using locate .."
if locate -l 1 -r "fftw3\.f03$" > /dev/null && locate -l 1 -r "libfftw3\(\.dll\)*\(\.a\|\.so\)$" > /dev/null ; then
   fftwpaths=`dirname $(locate -r "fftw3\.f03$") | uniq`
   echo ""
   echo "Which FFTW3 version should be used?"
   select d in $(echo $fftwpaths "other"); do
      if [ -n "$d" ]; then
         echo "$d selected"
         fftwpath=$d
         break
      fi
   done

   fftwlibs=`dirname $(locate -r "libfftw3\(\.dll\)*\(\.a\|\.so\)$") | uniq`
   echo ""
   echo "Which FFTW3 lib should be used?"
   select d in $(echo $fftwlibs "other"); do
      if [ -n "$d" ]; then
         echo "$d selected"
         fftwlib=$d
         break
      fi
   done
fi

if [ "$fftwlib" = "other" ] || [ "$fftwpath" = "other" ]; then
   echo ""
   echo "Automatic detection of FFTW3 failed."
   echo "You can either install it locally, or provide the path to the FFTW3 libary."
   read -e -p "Install under ~/.fftw? (y/n) " -i "n" dofftw
   case ${dofftw:0:1} in
      y|Y )
         FILE="fftw-3.3.4.tar.gz"
         dnfftw=""
         if [ ! -f $FILE ];
         then
            read -e -p "$FILE not found. Download from fftw.org? (requires wget) " -i "n" dnfftw
            case ${dnfftw:0:1} in
               y|Y )
                  command -v wget >/dev/null 2>&1 || { echo >&2 "I require wget but it's not installed.  Aborting."; exit 1; }
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
         fftwpath="$HOME/.fftw/include"
         fftwlib=`ls -d $HOME/.fftw/lib*`
         ;;
      * )
         echo "Please provide the path to the directory of the fftw3.f03 file:"
         read -e -i "/" fftwpath

         echo "Please provide the path to the directory of the FFTW3 libary (libfftw3):"
         read -e -i "/" fftwlib
   esac
fi
echo "FFTW_PATH = $fftwpath" > Src/make.fftwpath
echo "FFTWLIB = $fftwlib" >> Src/make.fftwpath
echo "yes"

echo -n "Checking ~/bin ..."
if [[ ! ":$PATH:" == *":$HOME/bin:"* ]]; then
   echo ""
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
      touch $conffile
      ;;
esac
