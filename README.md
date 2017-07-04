Molsim
==============

MOLSIM is an  integrated MC/molecular dynamics/Brownian dynamics simulation package created by [Per Linse in Lund](http://www.polyelectrolytes2002.fkem1.lu.se/). Version 6.4.7. of Molsim is beeing further developed and extended using this git.

Installation of Molsim
----------------------
*The following describes the installation of Molsim on the [RWTH-Compute Cluster](https://doc.itc.rwth-aachen.de/display/CC/Home).*

First Molsim clone the current repository. This requires having [set up your ssh key](https://git.rwth-aachen.de/help/ssh/README) at [git.rwth-aachen.de](https://git.rwth-aachen.de).
```shell
git clone git@git.rwth-aachen.de:pascal.hebbeker/Molsim.git
```
Now run the configure script. This will check some dependencies. Molsim requires FFTW 3.3.4. In can be install automatically within the configure script (Note: This might take some time). The configure script will also ask you for a version name. This version name will be appended to the executables of molsim (`molsim_ser.ver` instead of `molsim_ser`). Leave it blank for no special version name.
```shell
cd Molsim
./configure.sh
```
Now go to the Src directory and make Molsim:
```shell
cd Molsim/Src
make all mode=normal
```
For more details confer the [MOLSIM manual](http://www.polyelectrolytes2002.fkem1.lu.se/Molsim/Molsim.htm).

Usage
-----
Confer the [MOLSIM manual](http://www.polyelectrolytes2002.fkem1.lu.se/Molsim/Molsim.htm) for the standard commands. Confer the [wiki](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/home) for additional changes and features.

Additional Notes
----------------
* Cygwin
  * To install `FFTW` in cygwin you can either install it with the configure script (slow), or install the package `libfftw3-devel`. Either way you have to run `updatedb` afterwards to update the database of your file system.
