Appendix
========

# Namelist {#Namelist}
– A namelist consists of (i) one start statement, (ii), one end statement and (iii) a list (possibly
empty) between the start and end statements.

– The start statement consists of an ampersand (&) followed by the name of the namelist.

– The end statement consists of a slash (/).

– The list is composed of variables and their values separated with either by commas.

– The order of the variables and their assignments is normally irrelevant.

– However, if a variable is assigned values more than once, the last one takes place.

– Arrays are allocated either element by element, consecutively by listing values separated by
commas, or by a combination of these to ways.

– A value, say 5, appearing r times may be expressed as r*5.

– Note the order Fortran stores arrays elements: the left most index is running fastest and the
right most slowest.

– Comments are allowed and should proceed by an exclamation sign (!).

– Generally, if a program reads several namelists, they have to occur in the inputfile in the
same order as they are read from the program. Here, the order of the namelist is arbitrary.

– A namelist with no list (empty namelist) has to be specified.

– Namelists that are not read by the program do no harm.

Example: Consider the following declarations
```
character(4) :: title
integer(4) :: m, n
real(8) :: arr1(1:3), arr2(1:2,1:2), arr3(1:2,1:2)
```
The namelist  `nmlSystem`  below illustrates an element-by-element array assignment
```
&nmlSystem
title = 'head',
n = 10, m = 10,
n = 20,
arr1(1) = 1.5, arr1(2) = 2.5, arr1(3) = 2.5,
arr2(1,1) = 1.0, arr2(2,1) = 2.0, arr2(1,2) = 3.0, arr2(2,2) = 4.0,
arr3(1,1) = 1.0, arr3(2,1) = 2.0, arr3(1,2) = 3.0, arr3(2,2) = 4.0,
/
```
which could, for example, be shortened to
```
&nmlSystem
title = 'head',
n = 20,
m = 10,
arr1 = 1.5, 2*2.5,
arr2 = 1.0, 2.0, 3.0, 4.0,
arr3(1:2,1) = 1.0, 2.0, arr2(1:2,2) = 3.0, 4.0,
/
```

# Examples of namelists describing different objects {#ExampleNamelists}
Examples of namelist  `&nmlParticle` describing different types of objects are given in this appendix.

## Hard-sphere particles

`&nmlParticle`
> beginning of the namelist
`npt = 1,`
> number of particle types
``txpt = 'hs',``
> particle's label
`nppt = 100,`
> number of particles of this type
`natpt = 1,`
> number of atom types which this particle is composed of
``txat = 'hs',``
> atom's label
`massat= 1.0,`
> molar mass of the atom (in grams per mole)
`radat = 1.0,`
> hard-sphere radius of the atom
`naatpt(1,1) = 1,`
> number of 'hs' atoms in the 'hs' particle
``txaat(1,1) = 'hs',``
> label of 'hs' atom in the 'hs' particle
`/`
> end of the namelist

## Lennard-Jones particles

``&nmlParticle``
> beginning of the namelist
``npt = 1,``
> number of particle types
``txpt = 'lj',``
> particle's label
``nppt = 500,``
> number of particles of this type
``natpt = 1,``
> number of atom types which this particle is composed of
``txat = 'center',``
> atom's label
``massat= 40.0,``
> molar mass of the atom
``radat = 1.0,``
> hard-sphere radius of the atom
``sigat = 3.405,``
> Lennard-Jones parameter 'sigma'
``epsat = 0.99606555,``
> Lennard-Jones parameter 'epsilon'
``naatpt(1,1) = 1,``
> number of 'center' atoms in the 'lj' particle
``txaat(1,1) = 'lj',``
> label of 'center' atom in the 'lj' particle
``/``
> end of the namelist

## Lennard-Jones particles with an embedded dipole

``&nmlParticle``
> beginning of the namelist
``txelec= 'dip',``
> enable atoms to possess a static dipole
``npt = 1,``
> number of particle types
``txpt = 'LJ',``
> particle's label
``nppt = 100,``
> number of particles of this type
``natpt = 1,``
> number of atom types which this particle is composed of
``txat = 'LJ',``
> atom's label
``massat= 18.0,``
> molar mass of the atom
``radat = 0.5,``
> hard-sphere radius of the atom
``zat = 0.0,``
> atom's valency
``sigat = 2.88630,``
> Lennard-Jones parameter 'sigma'
``epsat = 1.97023,``
> Lennard-Jones parameter 'epsilon'
``naatpt(1,1) = 1,``
> number of 'LJ' atoms in the 'LJ' particle
``txaat(1,1) = 'LJ',``
> label of 'LJ' atom in the 'LJ' particle
``dipain(1,1,1) = 0.0000,``
> x, y, and z components of the static dipole of the first atom type of the first particle type
``0.0000, 0.10584,``
> end of the namelist
``/``

Note: Here, `dipain(1,1,1)` is a simplified version of `dipain(1:3,1,1)`.

## Nemo water (3-site water model with dipoles and polarizabilities)

``&nmlParticle``
> beginning of the namelist
``txelec= 'pol',``
> atoms possess charges, static and induced dipoles
``npt = 1,``
> number of particle types
``txpt = 'water',``
> particle's label
``nppt = 2,``
> number of particles of this type
``natpt = 2,``
> number of atom types which 'water' particles are composed of
``txat = 'o ','h ',``
> atom's labels; 1st atom type is oxygen
``massat= 16.0, 1.0,``
> molar mass of the atom
``radat = 0.0, 0.0,``
> hard-sphere radius of both atom types
``zat = -0.80100, 0.400500,``
> valency of oxygen and hydrogen atom
``naatpt(1,1) = 1, 2,``
> number of oxygen and hydrogen atom in water
``txaat(1,1) = 'o ','h ','h ',``
>  molecule atom's label
``rain(1,1,1) = 0.0, 0.0, -0.0656,``
> x, y, and z coordinates of oxygen
``rain(1,2,1) = 0.7572, 0.0, 0.5205,``
> x, y, and z coordinates of 1st hydrogen
``rain(1,3,1) = -0.7572, 0.0, 0.5205,``
> x, y, and z coordinates of 2nd hydrogen
``dipain(1,1,1) = 0.0000, 0.0000, -0.1299,``
> x, y, and z components of the static dipole of the oxygen atom
``dipain(1,2,1) = 0.0784, 0.0000, 0.0422,``
> x, y, and z components of the static dipole of the 1st hydrogen atom
``dipain(1,3,1) = -0.0784, 0.0000, 0.0422,``
> x, y, and z components of the static dipole of the 2nd hydrogen atom
``polain(1,1,1) = 0.6715, 0.6133, 0.7002, 0.0000, 0.0000, 0.0000,``
> xx, yy, zz, xy, xz, and zz components of the symmetric polarizability of the oxygen atom
``polain(1,2,1) = 0.2199, 0.0756, 0.1441, 0.0000, 0.1005, 0.0000,``
> xx, yy, zz, xy, xz, and zz components of the symmetric polarizability of the 1st hydrogen atom
``polain(1,3,1) = 0.2199, 0.0756, 0.1441, 0.0000, -0.1005, 0.0000,``
> xx, yy, zz, xy, xz, and zz components of the symmetric polarizability of the 2nd hydrogen atom
``/``
> end of the namelist

## Simple 1:1 elecrolyte

``&nmlParticle``
> beginning of the namelist
``npt = 2,``
> number of particle types
``txpt = 'ion1', 'ion2',``
> labels of both particle types
``nppt = 108, 108,``
> number of both particle types
``natpt = 1, 1,``
> number of atom types which particles are composed of
``txat = 'site1','site2',``
> labels of atom types
``massat= 23.0, 35.4,``
> molar mass of the atom (in grams per mole)
``radat = 2.0, 2.0,``
> hard-sphere radius of both atom types
``zat = 1.0, -1.0,``
> valency of both atom types
``naatpt(1,1) = 1,``
> number of 'site1' atoms in 'ion1' particle
``txaat(1,1) = 'site1',``
> label of 'site1' atom in 'ion1' particle
``naatpt(1,2) = 1,``
> number of 'site2' atoms in 'ion2' particle
``txaat(1,2) = 'site2',``
> label of 'site2' atom in 'ion2' particle
``/``
> end of the namelist

## Macroion + counterions

``&nmlParticle``
> beginning of the namelist
``npt = 2,``
> number of particle types
``txpt = 'macroion','ion',``
> labels of both particle types
``nppt = 1, 60,``
> number of particles of each particle types
``natpt = 1, 1,``
> number of atom types which particles are composed of
``txat = 'mic', 'ion',``
> labels of atom types
``massat= 460.0, 23.0,``
> molar mass of the atom (in grams per mole)
``zat = -60, 1,``
> valency of both atom types
``radat = 20.0, 2.0,``
> hard-sphere radius of both atom types
``naatpt(1,1) = 1,``
> number of 'mic' atoms in 'macroion' particle
``txaat(1,1) = 'macroion',``
> label of 'macroion' atom in 'macroion' particle
``naatpt(1,2) = 1,``
> number of 'ion' atoms in 'ion' particle
``txaat(1,2) = 'ion',``
> label of 'ion' atom in 'ion' particle
``/``
> end of the namelist

##  Polyions + counterions

``&nmlParticle``
> beginning of the namelist
``nct = 1,``
> number of chain types
``txct ='100-mer',``
> label of the chain type
``ncct = 10,``
> number of chains in the system
``npptct(1,1) = 100, 0``
> number of particles in the '100-mer' is equal to 100 of type 'pe' and zero of type ion'
``npt = 2,``
> number of particle types
``txpt = 'pe', 'ion',``
> labels of both particle types
``nppt = 1000,1000,``
> number of particles of each particle types
``natpt = 1, 1``
> number of atom types which particles are composed of
``txat = 'bead', 'ion',``
> labels of atom types
``massat= 10.0,10,``
> molar mass of the atom (in grams per mole)
``radat = 2.0, 2.0,``
> hard-sphere radius of both atom types
``zat = 1.0, -1.0,``
> valency of both atom types
``naatpt(1,1) = 1``
> number of 'bead' atoms in 'pe' particle
``txaat(1,1) = 'bead',``
> label of 'bead' atom in 'pe' particle
``naatpt(1,2) = 1,``
> number of 'ion' atoms in 'ion' particle
``txaat(1,2) = 'ion',``
> label of 'ion' atom in 'ion' particle
``/``
> end of the namelist

Note: Here `npptct(1,1)` is a simplified version of `npptct(1:2,1)`.

## Diamond-like polyelectrolyte gel + counterions

``&nmlParticle``
> beginning of the namelist
``lclink=.true.,``
> flag to enable crosslinking
``maxnbondcl= 4, 2, 0,``
> maximum number of crosslinks for types of particles
``nct = 1,``
> number of chain types
``txct ='strand'``
> label of chain type
``ncct = 16,``
> number of chains of type 'strand'
``npptct(2,1) = 10,``
> number of particles of type 'strand' in the chain type 'strand'
``npt = 3,``
> number of particle types
``txpt='node','strand','countion',``
> labels of all particle types
``nppt = 8, 160, 168,``
> number of particles of each particle types
``natpt = 1, 1, 1``
> number of atom types which particles are composed of
``txat='bead1','bead2','countion'``
> labels of atom types
``radat = 2.0, 2.0, 2.0,``
> hard-sphere radius of both atom types
``zat = 1.0, 1.0, -1.0,``
> valency of all three atom types
``naatpt(1,1) = 1,``
> number of 'bead' atoms in 'node' particle
``txaat(1,1) = 'bead',``
> label of 'bead' atom in 'node' particle
``naatpt(1,2) = 1,``
> number of 'bead2' atoms in 'strand' particle
``txaat(1,2) = 'bead2',``
> label of 'bead2' atom in 'strand' particle
``naatpt(1,3) = 1,``
> number of 'countion' atoms in 'countion' particle
``txaat(1,3) = 'countion',``
> label of 'countion' atom in 'countion' particle
``/``
> end of the namelist

Note: Here `npptct(1,1)` is a simplified version of `npptct(1:2,1)`.
## Bottle-brushes polymers with uneven side-chain distribution

``&nmlParticle``
> beginning of the namelist
``lclink =.true.,``
> flag to enable crosslinking
``maxnbondcl = 1, 1,``
> maximum number of crosslinks for types of particles
``ngen = 1,``
> number of generations of hierarchichal structure
``ictgen(0) = 1,``
> chain type of 0th generation (main chain)
``ictgen(1) = 2,``
> chain type of 1st generation (side chain)
``nbranch = 5,``
> number of branches on the main chain
``ibranchpbeg = 1,``
> main chain particle number of the first branch point
``ibranchpinc = 1,``
> segment increment for branch points on the main chain
``nct = 2,``
> number of chain types
``txct ='10-mer1','5-mer2',``
> label of each chain type
``ncct = 2, 10,``
> number of chains of each type
``npptct(1,1) = 10, 0,``
> number of particles of each type in the 1st chain
``npptct(1,2) = 0, 5,``
> number of particles of each type in the 2nd chain
``npt = 2,``
> number of particle types
``txpt = 'bead1', 'bead2',``
> labels of both particle types
``nppt = 20, 50,``
> number of particles of each particle types
``natpt = 1, 1,``
> number of atom types which particles are composed of
``txat = 'site1', 'site2',``
> labels of atom types
``massat= 10.0, 10.0,``
> molar mass of the atom (in grams per mole)
``radat = 2.0, 2.0,``
> hard-sphere radius of both atom types
``naatpt(1,1) = 1,``
> number of 'site1' atoms in 'bead1' particle
``txaat(1,1) = 'site1',``
> label of 'site1' atoms in 'bead1' particle
``naatpt(1,2) = 1,``
> number of 'site2' atoms in 'bead2' particle
``txaat(1,2) = 'site2',``
> label of 'site2' atoms in 'bead2' particle
``/``
> end of the namelist

Note: Here `npptct(1,1)` is a simplified version of `npptct(1:2,1)` etc.
```
particle id
  1        5        10  11      15         20
  x x x x x x x x x x    x x x x x x x x x x
 21  O  O  26  O  41     46  O  O  51  O  66
   O O      O            O O    O
   O O ...  O            O O .. O
   O O      O            O O    O
   O O      O            O O    O
```

## Third generation dendrimers

``&nmlParticle``
> beginning of the namelist
``lclink =.true.,``
> flag to enable crosslinking
``maxnbondcl = 2, 2, 2, 1,``
> maximum number of crosslinks for types of particles
``ngen = 3``
> number of generations of hierarchichal structure
``ictgen(0) = 1,``
> chain type of 0th generation
``ictgen(1) = 2, 3, 4,``
> chain type of 1st, 2nd, and 3rd generation
``nbranch = 2, 2, 2,``
> number of branches of 1st, 2nd and 3rd generation
``ibranchpbeg = 1, 4, 4,``
> segment number of the first branch point of 1st, 2nd, and 3rd generation
``ibranchpinc = 0, 0, 0,``
> segment increment for branch points of 1st, 2nd, and 3rd generation
``nct = 4,``
> number of chain types
``txct ='1-mer1', '4-mer2', '4-mer3', '4-mer4'``
> label of each chain type
``ncct = 1, 2, 4, 8,``
> number of chains of each type
``npptct(1,1) = 1, 0, 0, 0,``
> number of particles of type 'bead1' in the chain type '1-mer1'. Chain '1-mer1' is composed only of 'bead1' particles.
``npptct(1,2) = 0, 4, 0, 0,``
> number of particles of type 'bead2' in the chain type '4-mer2'. Chain '4-mer2' is composed only of 'bead2' particles.
``npptct(1,3) = 0, 0, 4, 0,``
> number of particles of type 'bead3' in the chain type '4-mer3'. Chain '4-mer3' is composed only of 'bead3' particles.
``npptct(1,4) = 0, 0, 0, 4,``
> number of particles of type 'bead4' in the chain type '4-mer4'. Chain '4-mer4' is composed only of 'bead4' particles.
``npt = 4,``
> number of particle types
``txpt = 'bead1', 'bead2', 'bead3', 'bead4',``
> labels of all particle types
``nppt = 1, 8, 16, 32,``
> number of particles of each particle types
``natpt = 1, 1, 1,``
> number of atom types which particles are composed of
``txat = 'site1', 'site2', 'site3', 'site4'``
> labels of atom types
``massat= 10., 10., 10., 10.,``
> molar mass of the atom (in grams per mole)
``radat = 2.0, 2.0, 2.0, 2.0,``
> hard-sphere radius of both atom types
``naatpt(1,1) = 1,``
> number of 'site1' atoms in 'bead1' particle
``txaat(1,1) = 'site1',``
> label of 'site1' atoms in 'site1' particle
``naatpt(1,2) = 1,``
> number of 'site2' atoms in 'bead2' particle
``txaat(1,2) = 'site2',``
> label of 'site2' atoms in 'site2' particle
``naatpt(1,3) = 1,``
> number of 'site3' atoms in 'bead3' particle
``txaat(1,3) = 'site3',``
> label of 'site3' atoms in 'site3' particle
``naatpt(1,4) = 1,``
> number of 'site4' atoms in 'bead4' particle
``txaat(1,4) = 'site4',``
> label of 'site4' atoms in 'site4' particle
``itestpart = 10,``
> test output of chain pointers
``/``
> end of the namelist

Note: Here `npptct(1,1)` is a simplified version of `npptct(1:2,1)` etc.
```
particle id
     30    33   49  46
   13+ # # # #   # # # # +21
  26# +             + #42
   #   +          +  #
 29# #  10+5  2 1 6  9+18  # #45
     O O O O x O O O O
 37# #  14+        +22 #  #53
    #  +         + #
   34# +          + #50
    17+ ####     #### +25
     38  41    57  54
```

# Suggestion or Problems? {#SuggestionOfProblems}

If you have suggestions or problems with Molsim, open an issue [online](https://github.com/joakimstenhammar/molsim/issues/). Please add the following information to help us help you:

* MOLSIM version number (check by running `molsim --version`)
* Host computer (name and configuration)
* Suggestion or problem in as much detail as possible
* Any Input/Cnf files needed to reproduce the error
