User-provided routines
======================
The file  moluser.F90 contains routines that are conditionally called and which might be modified by
the user. This allows the software to be extended without having to modify the core of it.

The existing code in the distributed file  moluser.F90 reflects the current use of it. The procedure to
insert a new routine is as follows:

 –  At the appropriate position in driver moluser, extend the test of  txUser and write the
appropriate call to the new subroutine

 – Add your new subroutine in  moluser.F90

 –  Compile

The file  moluser.F90 contains six drivers and several routines specifying different types of tasks.
The six drivers are  PotentialUser ,  SetParticleUser ,  DumpUser ,  GroupUser ,  StaticUser , and  ImageUser .

# PotentialUser {#PotentialUser}
Routine PotentialUser is the driver of user-provided routines for two-body potentials. It calls
routines specifying various potentials on the basis of the value of txpot(iptjpt), which is specified in
namelist nmlPotential. The procedure to install a new user-provided potential is as follows:

 – Extend the test of  txpot(iptjpt) , and write the appropriate call to a new routine in routine
PotentialUser.

 – Add your new routine to file  moluser.F90.

# SetParticleUser {#SetParticleUser}
Routine  SetParticleUser is the driver of user-provided routines for generating a start configuration.
It calls routines specifying types of initial particle positions on the basis of the value of
txsetconf(ipt) , which is specified in namelist nmlSetConfiguration. The procedure to install a new
user-provided start configuration is as follows:

 –  Extend the test of  txsetconf(ipt) and write the appropriate call to a new routine in the
routine SetParticleUser.

 – Add your new routine to file  moluser.F90.

# DumpUser {#DumpUser}
Routine DumpUser is the driver of user-provided routines for dumping data. It calls routines
specifying different types of dumping. Generally, all these user-provided routines are called if
variable  ldumpuser , specified in namelist  nmlDump , is true. The procedure to install a new dump
analysis is as follows:

 –  Write the appropriate call to the new routine in the routine  DumpUser.

 – Add your new routine to file  moluser.F90

# GroupUser {#GroupUser}
Routine  GroupUser is the driver of user-provided routines for particle group division. It calls
routines specifying different division of particles into groups on the basis of the values of variables
ref or field [later stored in  txappl(m) (m = 1 for ref and m = 2 for field particles)], which are
specified in namelist  nmlGroup . The procedure to install a new group division procedure is as
follows:

 –  Extend the test of  txappl(m) (m = 1 for ref and m = 2 for field particles) and write the
appropriate call to the new routine in driver  GroupUser.

 –  Add your new routine to the file  moluser.F90.

Follow closely the pattern of one of the other user-provided routines for particle group division.
After the statement 'if (iStage = iWriteInput) then', variables  ngr ,  igrpt , and  txgrr have to be
specified.  ngr denotes the number of groups,  igrpt the order number of the particle type of particles
that might belong to this group, and  txgrr a text label of the group. After 'if (iStage = iSimStep)
then', the variable  igr (within the particle loop) has to be assigned the order number of the group to
which the particle should belong (zero if the particle does not belong to any group). The last two
lines of the particle loop should be the same as in the other cases.

# StaticUser {#StaticUser}
Routine  StaticUser is the driver of user-provided routines for static analysis. It calls routines
specifying different types of analysis. Generally, all these user-provided routines are called if
variable lstaticuser, specified in namelist  nmlStatic , is true. The procedure to install a new static
analysis routine is as follows:

 –  Write the appropriate call to the new routine in the routine  StaticUser.

 – Add your new routine to the file  moluser.F90.

# ImageUser {#ImageUser}
Routine  ImageUser is the driver of user-provided image data writing routines. It calls routines
specifying different types of such routines. Generally, all these user-provided routines are called if
variable limageuser, specified in namelist  nmlImage , is true. The procedure to install a new image
data writing routine is as follows:

 –  Write the appropriate call to the new routine in the routine  ImageUser.

 – Add your new routine to the file  moluser.F90.
