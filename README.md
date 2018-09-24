Molsim
==============

MOLSIM is an  integrated MC/molecular dynamics/Brownian dynamics simulation package created by [Per Linse in Lund](http://www.polyelectrolytes2002.fkem1.lu.se/). Version 6.4.7. of Molsim is beeing further developed and extended using this git.

Obtaining the Code
------------------
There are two possibilites of how to obtain the code. You can either simply download the tarball of the code, or clone the whole repository.
### Downloading the tarball
Download the tarball from [here](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/repository/archive.tar.gz?ref=master) and save it in the directory of your choice. Afterwards, navigate to that directory and exctract the source code with
```
tar -zxf <name of the tar file>
```
You might want to rename the directory to some more resonable name.

### Clone the Molsim repository

This requires having [set up your ssh key](https://git.rwth-aachen.de/help/ssh/README) at [git.rwth-aachen.de](https://git.rwth-aachen.de). Simply run
```shell
git clone git@git.rwth-aachen.de:pascal.hebbeker/Molsim.git
```
Installation of Molsim
----------------------

Navigate into the Molsim directroy and run the configure script. This will check some dependencies. Molsim requires FFTW 3.3.4. In can be install automatically within the configure script (Note: This might take some time). The configure script will also ask you for a version name. This version name will be appended to the executables of molsim (`molsim_ser.ver` instead of `molsim_ser`). Leave it blank for no special version name.
```shell
cd Molsim
./configure.sh
```
The configure script tries to locate the FFTW libary. If you want to customize the path of the libary, modify the `Src/make.fftwpath` file. Additionally it will select which compiler to use. To costumize the compiler which is to be used change the `Src/make.arch` file.

Now make Molsim:
```shell
make all
```
For more details confer the [MOLSIM manual](documentation.pdf).

Usage
-----
Confer the [MOLSIM manual](documentation.pdf) for the standard commands.

Additional Notes
----------------
* Cygwin
  * To install `FFTW` in cygwin you can either install it with the configure script (slow), or install the package `libfftw3-devel`. Either way you have to run `updatedb` afterwards to update the database of your file system.
