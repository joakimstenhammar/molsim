External units
==============
A file is specified by its name and type, which are separated by a dot. The file types are used in a
predetermined way to describe the use of the file and leaves file name to label the project. Within a
project the file name is independent of program.

FPOS, FORI, FLIV, FANV, FFOR, FTOR, and FIDM are collectively called dump files.

<table>
<tr>
<th>File type</th>
<th>Generic name</th>
<th>Content</th>
<th>binary</th>
</tr>
<tr>
<td>
in<br>
out<br>
cnf<br>
lib<br>
list<br>
img<br>
user<br>
group<br>
pos<br>
ori<br>
liv<br>
anv<br>
for<br>
tor<br>
idm<br>
</td>
<td>
FIN<br>
FOUT<br>
FCNF<br>
FLIB<br>
FLIST<br>
FIMG<br>
FUSER<br>
FGROUP<br>
FPOS<br>
FORI<br>
FLIV<br>
FANV<br>
FFOR<br>
FTOR<br>
FIDM<br>
</td>
<td>
Input data<br>
Output data<br>
Configuration and average data<br>
Potential library<br>
List data<br>
Image data<br>
At the disposal to the user<br>
Particle group data<br>
Box lengths and positions of particles<br>
Orientations of particles<br>
Linear velocities (box frame) of particles<br>
Angular velocities (principal frame) of particles<br>
Forces (box frame) of particles<br>
Torques (box frame) of particles<br>
Induced dipole moments (box frame) of particles<br>
</td>
<td>
   <br>
   <br>
binary<br>
   <br>
   <br>
   <br>
   <br>
   <br>
binary<br>
binary<br>
binary<br>
binary<br>
binary<br>
binary<br>
binary<br>
</td>
</tr>
</table>

\image latex molsim_in_out.png
**Flow chart of the use of external files.** Bold arrows denote files that are always required or
generated, while thin arrows denote files that might be required or generated, depending on input
variables.

