Additional Descriptions
=======================

# Networks {#Networks}

## Overview
### Non-periodic network structures in Molsim
The aim of this branch  is the "proper" implementation of polymer networks as an integral part of the simulation software Molsim. These networks are called non-periodic networks in order to distinguish them from their macroscopic counter parts (i.e. macrogels). This follows Molsim's philosophy of having different objects in a hierarchy of different levels.

| Object type             | Hierarchy Level | Object type                     | Hierarchy Level |
| ---                     | ---             | ---                             | ---             |
| Atoms                   | 1               |                                 |                 |
| Particles               | 2               |                                 |                 |
| Chains                  | 3               |                                 |                 |
| Hierarchical Structures | 4               | Non-periodic Network Structures | 4               |

Within this scheme, atoms are part of particles, which in turn may be connected to form a chain. These can be interconnected by cross-links to form more complex structures, as e.g. hierarchical structures or network-like structures. Note, that for all object types on the left hand side of the table, a programming infrastructure in form of various characteristic parameters and (pseudo-)pointers is established already, thus, leading to a high programming flexibility, whenever one wishes to expand the functionality of the software with respect to one of these object types. On the other hand, non-periodic network structures have been introduced to Molsim later and therefore a programming infrastructure is not given. Besides the mentioned characteristic parameters and (pseudo-)pointers, a series of subroutines in Molsim handles the statistical treatment of the left hand side object types. Within the scope of this branch the corresponding parameters and functionalities shall be established.

### Model description
The aim of the implemented network model is to simulate the properties of non-periodic polymer networks of spherical shape (i.e. microgels). The model network comprises cross-linking particles (so called nodes) and chain particles. The nodes are set on the carbon positions of the diamond lattice (space group Fd3m, no. 227). The cubic lattice is cropped to obtain a spherical shape by defining a cutoff radius - all nodes within the cutoff radius will be kept; all nodes outside of the cutoff radius will be discarded. The remaining nodes are then connected via polymer chains, which consist of a number of chain particles per chain. Within the model the so called dangling chains are included. They are connected to only one cross-link and thus they form the outer periphery of the modelled gel. Via the input parameter `nnwt` it is possible to define the number of network types present in the system. Different network types may for example differ in their number (`nnwnwt(inwt)`) and size (`rnwt(inwt)`). Further, different types of particles may be used as building blocks and hence, polymer networks with different physical properties may be obtained. Within this scheme it is important to note, that always one particle type forms the nodes of one certain network type. Analogues, one chain type forms the chains of a certain network type. This is controlled via the parameters `iptclnwt(inwt)` (i.e. the particle type forming the cross-links of network type `inwt`) and `ncctnwt(ict,inwt)` (number of chains of a chain type `ict` belonging to one network of network type `inwt`). The size of the networks cannot be chosen in a continous manner, as the nodes are set on discrete diamond lattice positions, hence the generated networks are available with discrete numbers of nodes.

## Details of implementation
### Network Configuration
### General
The possibility to set finite networks has already been possible prior to the here undertaken implementations by means of the subroutine `SetNetwork`. Within the branch Networks the already existing routine `SetNetwork` has been adapted to the new nomenclature regarding network related parameters (see below). Both, in the subroutine `SetNetwork` and in the general output about the system parameters, the output about networks has been refined.
### Shift of cropping sphere
The parameter `shiftnwt(1:3,1:nnwt)` was introduced. It enables one to obtain topology-modified networks, which mainly differ in how the structure (*i.e.* dangling chains and interconnectivity of the outer layer) of their periphery looks like. In the process of generating the network structure a cubic diamond lattice will be span and cropped by a spherical cut-off. The default network will be set having the (0,0,0) position of the diamond lattice unit cell in the center of the cropping sphere. This corresponds to `shiftnwt(1:3,inwt) = 0.0` of a certain network type `inwt`. By setting 'shiftnwt(1:3,inwt)' the center of the cropping sphere may be shifted to any point of the unit cell. Note, that `shiftnwt` accepts its values in the range of 0.0 to 1.0, respectively. Higher values will work aswell, but due to the underlying periodic boundary conditions which the diamond lattice is subject to, `shiftnwt(1:3,inwt) = 0.0` is equivalent to `shiftnwt(1:3,inwt) = 1.0` and so on.

### Network Properties
Similar to the possibility to evaluate characteristic properties of chain objects in Molsim, a new subroutine has been implemented, namely `CalcNetworkProperty`. Within the subroutine `CalcNetworkProperty` the following properties are evaluated:
* `ro(3)`: center of mass
* `rg2`: radius of gyration squared
* `rg2x`: radius of gyration squared projection on x-axis
* `rg2y`: radius of gyration squared projection on y-axis
* `rg2z`: radius of gyration squared projection on z-axis
* `rg2s`: square extention along principal axes (smallest)
* `rg2m`: square extention along principal axes (middle)
* `rg2l`: square extention along principal axes (largest)
* `asph`: asphericity (JCP 100, 636 (1994))
* `theta(3)`: angle between x-, y- and z-axis, respectively, and the axis of the network's largest extension
* `alpha`: degree of ionization (for titrating systems)
* `eivr(3,3)`: normalized eigenvectors of the principal frame

### Averages of Network Properties
Within the subroutine `NetworkAver` the averages and precisions with regards to the properties available by `CalcNetworkProperties` are formed (except for `ro(3)` and `eivr(3,3)`). The obtained network quantities will be written to the Output-File when networks are present in the system.
Example output:

```
*************************************************************************************
*                                    network quantities                             *
*************************************************************************************
quantity                         average  precision  fluctuation  precision  stat eff
--------                         -------  ---------  -----------  ---------  --------
<r(g)**2>**0.5               =  75.32408   0.00915       1.77217    0.06023        2.
<r(g)**2_x>**0.5             =  43.51678   0.00452       1.36440    0.09226        0.
<r(g)**2_y>**0.5             =  43.49056   0.00458       1.37475    0.14198        0.
<r(g)**2_z>**0.5             =  43.45776   0.01588       1.66017    0.01477        3.
smallest rms mom. p.a.       =  43.39747   0.01547       1.77924    0.10399        2.
intermediate rms mom. p.a.   =  43.47212   0.01324       1.65132    0.11078        2.
largest rms mom. p.a.        =  43.59531   0.01354       1.74615    0.07283        1.
<asphericity>                =   0.00001   0.00000       0.00000    0.00000        1.
<xtheta>                     =  72.47060   8.72678      44.88033    3.26752        0.
<ytheta>                     = 101.41857   3.97221      38.12409    4.66430        0.
<ztheta>                     =  78.86201   2.54458      19.66774    2.88301        0.

asphericity (<2:nd moments>) =   0.00001
```
Here, the asphericity is being calculated according to [Zifferer, G. and Olaj, O. F., Journal of Chemical Physics 1994, 100, 636](http://dx.doi.org/10.1063/1.466926). The second displayed asphericity (`asphericity (<2:nd moments>)`) is being directly calculated from the respective average of the principal moments of a macrostep or the whole simulation.

### Network Distribution Functions
The functionality of finite network-related simulations has been extended by allowing for statistical analysis of network properties by means of distribution functions.

###  Network Property Distribution Functions
The following network property distribution functions are available by means of the setting in the `nmlNetworkDF`.

| type | label      | description                                              |
| ---  | ---        | ---                                                      |
| 1    | `rg`       | network radius of gyration distribution                  |
| 2    | `asph`     | network asphericity distribution                         |
| 3    | `alpha`    | network degree of ionization distribution                |
| 4    | `rgchain`  | chain radius of gyration distribution of network chains  |
| 5    | `reechain` | chain end-to-end distance distribution of network chains |

The described type `i` refers to the array index to be used in the respective `vtype(i)` construct.

### Radial Network Distribution Functions
The following radial network distribution functions are available by means of the setting in the `nmlNetworkRadialDF`.

| type | label     | description                          |
| ---  | ---       | ---                                  |
| 1    | `rpart`   | radial particle number distribution  |
| 2    | `rdens`   | radial particle density distribution |
| 3    | `rgchain` | radial chain radius of gyration      |
| 4    | `q`       | radial sum of all charges            |
| 5    | `qcum`    | radial cumulated sum of all charges  |
| 6    | `alpha`   | radial degree of ionization          |
| 7    | `rchain`  | radial chain number distribution     |

The described type `i` refers to the array index to be used in the respective `vtype(i)` construct.

### Network Generation Groups
A new way of forming groups within networks has been established. Within this grouping scheme the chains are numbered starting from the dangling chains at the periphery of the network and taking the shortest way to the core chains of the network structure. A schematic representation of the numbering is presented below.
```
           | 1
           o
           | 2
           o
 1   2   3 | 3   2   1
 - o - o - o - o - o -
           | 3
           o
           | 2
           o
           | 1
```
Note, that this scheme works only for structures set with `txsetconf = 'network'`. The cross-links of the network will be assigned to have no group. In order to obtain this grouping set either `ref` or `field` in the `nmlGroup` to `networkgenerations`.

## Testin
### fnw.nwt2.mc.in
 The testin-file was added to the in-directory. The corresponding out-file was stored in the Save-directory. The documents "todiff.txt" and "goall.sh" were expanded by the new input-file.
### fnw.nwgen.df.mc.in
 The testin-file was added to the in-directory. The corresponding out-file was stored in the Save-directory. The documents "todiff.txt" and "goall.sh" were expanded by the new input-file.

# Network variables
### Internal Variables
Besides of the external network related parameters (i.e. input variables) a number of internal variables has been introduced in order to set up the programming infrastucture of non-periodic networks. All variables are declared in `mol.F90`.

| Keyword                   | Type            | Description                                                                |
| ---                       | ---             | ---                                                                        |
| `nnw`                     | `integer`       | Total number of networks                                                   |
| `nctnwt(1:nnwt)`          | `integer`       | Number of different chain types in network type `inwt`                     |
| `ncnwt(1:nnwt)`           | `integer`       | Number of chains per network of network type `inwt`                        |
| `npnwt(1:nnwt)`           | `integer`       | Number of particles per network of network type `inwt`                     |
| `nclnwt(1:nnwt)`          | `integer`       | Number of cross-links per network of network type `inwt`                   |
| `nnwtnwt`                 | `integer`       | Number of different network type pairs                                     |
| `txnwtnwt(1:nnwtnwt)`     | `character(21)` | Name of network type pair `inwtnwt`                                        |
| `massnwt(1:nnwt)`         | `real`          | Mass per network of network type `inwt`                                    |
| `massinwt(1:nnwt)`        | `real`          | Inverse mass per network of network type `inwt`                            |
| `lnetwork`                | `logical`       | `.true.`if networks are used                                               |
| `lptnwt(1:mnpt,1:nnwt)`   | `logical`       | `.true.` if particle type `ipt` is part of networks of network type `inwt` |
| `lpnnwn(1:np,1:nnw)`      | `logical`       | `.true.` if particle `ip` is part of network `inw`                         |
| `npweakchargenwt(1:nnwt)` | `integer`       | Number of titratable charges per network of network type `inwt`            |

### Network pointer
| Keyword                       | Type      | Description                                                                                  |
| ---                           | ---       | ---                                                                                          |
| `inwtnwn(1:nnw)`              | `integer` | In: network number `inw`, Out: its network type `inwt`                                       |
| `inwtct(1:nct)`               | `integer` | In: chain type `ict`, Out: its network type `inwt`                                           |
| `inwtcn(1:nc)`                | `integer` | In: chain number `ic`, Out: its network type `inwt`                                          |
| `inwncn(1:nc)`                | `integer` | In: chain number `ic`, Out: its network number`inw`                                          |
| `inwnnwt(1:nnwt)`             | `integer` | In: network type `inwt`, Out: its first network `inw`                                        |
| `inwtnwt(1:nnwt,1:nnwt)`      | `integer` | In: two network types `inwt/jnwt`, Out: its network type pair number `inwtnwt`               |
| `icnclocnwn(1:ncnwt,1:nnw)`   | `integer` | In: local chain number `icloc` and newtwork number `inw`, Out: its chain number `ic`         |
| `ipncllocnwn(1:nclnwt,1:nnw)` | `integer` | In: local cross-link number `iclloc` and network number `inw`, Out: its particle number `ip` |
| `ipnplocnwn(1:npnwt,1:nnw)`   | `integer` | In: local particle number `iploc` and network number `inw`, Out: its particle number `ip`    |

# Complexation-Analysis {#ComplexationAnalysis}
## Complexation Analysis
Interparticle complexation can be analyzed using Molsim.  To use it set

```fortran
lstaticuser = .true.
```
in `nmlStatic`. Additionallly the first 4 characters of `txuser` in `nmlSystem` have to be `comp`. The control of the complex-statistics is given by the  namelist `nmlComplexation`:
```fortran
 &nmlComplexation
  rcut_complexation = ...,
  lComplexFraction = ....,
  lClusterDF = ...,
 /
```
control over the ComplexDF routine is given by the namelist `nmlComplexDF`:
```fortran
 &nmlComplexDF
   vtype = ...
 /
```
where `rcut_complexation` defines the maximum distance for two particles to form a complex; `lComplexFraction` (`logical`) switches on the calculation of the fraction of complexed particles; `lClusterDF` (`logical`) switches on the calculation of the size distribution of complexed formed by the complexation.

* Testfile: `complexation.mc.in`

### Fraction of Complexation
This routine return the fraction of complexed particles of each particle-particle commbination. It return the fraction of `X` beads which are closer than `rcut_complexation` to a `Y` bead as `w(cmplx): X - Y`. E.g.:

```fortran
quantity           average  precision  fluctuation  precision  stat eff
--------           -------  ---------  -----------  ---------   -------
w(cmplx): A - A    0.30000    0.00000      0.00000    0.00000        0.
w(cmplx): A - B    0.20000    0.00000      0.00000    0.00000        0.
w(cmplx): B - A    0.40000    0.00000      0.00000    0.00000        0.
w(cmplx): B - B    0.50000    0.00000      0.00000    0.00000        0.
```
means that 30% of the A beads are complexed to other A beads;  20% of the A beads are complexed to a `B` bead; 40% of the B beads form a complex with an A bead and 50% of the B beads form a Complex with another B bead.

### Complexation Distribution
This routine return the distribution of the fraction of complexed particles of each particle-particle commbination. It return the fraction of `X` beads which are closer than `rcut_complexation` to a `Y` bead. To control of the distribution function the `nmlComplexDist` is used. Additionally also the distribution of the complexation rate of the individual chains can be calculated.

### Segment Complexation
This routine return the fraction of complexation of each segment along the chains of the system (averaged over all chains of one type).

### Size distribution of the binary clusters
Here the size of clusters formed by neighbouring complexed particles of two different types are measured. A cluster is formed only by connections between particles of different particle type. E.g.
```
A   A
    |
    B
    |
    A
    |
A - B
```
forms one cluster of size 5  and one of size 1 (single uncomplexed A particle). In contrast
```
A - B
    |
    A
     
    A
    |
A - B
```
forms two clusters of size 3.

As an output the size distribution (average number of cluster of a specific size) are given for all combinations of two different particle types.

# Advanced Configurations {#Configurations}
The support for branched, random and repeating structures is improved:
## Branched Structures
* When setting a hierarchical structure using the subroutine `SetHierarchical` the inability to set one particle will not lead to an abortion of the whole program. Instead another it is attempted `ntrydef` (default 100) times to set the whole structure. If one particle cannot be set the whole structure is reset.
* Testfile: `hierarchical.mc.in`

# Copolymer Sequence {#CopolymerSequence}
## Overview
In general chains consisting of more than one particle type (i.e. copolymers) may be simulated using Molsim. Different ways to generate such chains are triggered by the parameter `txcopolymer`. Existing modes comprise alternating, and block copolymers. Within the scope of this branch the functionality to simulate copolymers with random, repeating or irregular or highly specific monomer distribution has been added (sequence).

## Usage
### Sequence
The input of the copolymer sequence has been realized by introducing the parameter `iptsegct(iseg,ict)`, which specifies the particle type `ipt` of chain segments `iseg` of chains of type `ict`.
First of all the parameter `txcopolymer(ict)` has to be set to `sequence` for the chain type in question in the namelist `nmlParticle`:
```
&nmlParticle
 ...
 txcopolymer(ict) = 'sequence' ,
 ...
/
```
where `ict` should be replaced by the corresponding chain type number. The parameter `iptseg(iseg,ict)` is then being set in the namelist `nmlCopolymerSequence`. For example
```
&nmlCopolymerSequence
 iptsegct(1,ict) = ipt , iptsegct(2,ict) = ipt , ... , iptsegct(npct,ict) = ipt ,
 ...
/
```
where again `ict` should be substituted by the corresponding chain type number, `ipt` should be replaced by the intended particle type number and `npct` corresponds to the number of particles per chain of chain type `ict`.


### Repeating Copolymers
Copolymers with a defined repeating block structure can be generated. It consists of repeating units, each consisting of blocks of one particle type. 
* in `nmlParticle` set the parameters of the chain type `ict` which is to be repeating

```fortran
txcopolymer(ict) = "repeating",
```
and set

``` fortran
nblockict(ict)
```
to the number of blocks within a repeating unit to the intended number. 
* The detailed repeating structure is given by `nmlRepeating`:

```fortran
&nmlRepeating
rep_iblock_ict( iblock , ict)%pt = ... ,
rep_iblock_ict(iblock , ict)%np = ... ,
/
```
where the particle type `%pt` and number of particles `%np` of each `iblock` for each`ict` is defined.

If not enough particles of a certain type are given for a chain to fulfill the regular structure, the repeating structure is continued, but only filling the blocks with the remaining particles.

##### Example:
One chain with the following input (excerpt)

```fortran
 &nmlParticle
  nct   = 1,
  ncct  = 1,
  npptct(1:3,1) = 4,5,10,
  nblockict(1) = 3,
  txcopolymer='repeating',
 /
 &nmlRepeating
 rep_iblock_ict(1:3,1)%pt = 2, 3, 1
 rep_iblock_ict(1:3,1)%np = 2, 3, 1
 /
```
will generate such a repeating structure:

<table>
  <tr>
    <th>Segment</th>
    <td>1</td>
    <td>2</td>
    <td>3</td>
    <td>4</td>
    <td>5</td>
    <td>6</td>
    <td>7</td>
    <td>8</td>
    <td>9</td>
    <td>10</td>
    <td>11</td>
    <td>12</td>
    <td>13</td>
    <td>14</td>
    <td>15</td>
    <td>16</td>
    <td>17</td>
    <td>18</td>
    <td>19</td>
  </tr>
  <tr>
    <th>Particle Number</th>
    <td>5</td>
    <td>6</td>
    <td>10</td>
    <td>11</td>
    <td>12</td>
    <td>1</td>
    <td>7</td>
    <td>8</td>
    <td>13</td>
    <td>14</td>
    <td>15</td>
    <td>2</td>
    <td>9</td>
    <td>16</td>
    <td>17</td>
    <td>18</td>
    <td>3</td>
    <td>19</td>
    <td>4</td>
  </tr>
  <tr>
    <th>Particle Type</th>
    <td>2</td>
    <td>2</td>
    <td>3</td>
    <td>3</td>
    <td>3</td>
    <td>1</td>
    <td>2</td>
    <td>2</td>
    <td>3</td>
    <td>3</td>
    <td>3</td>
    <td>1</td>
    <td>2</td>
    <td>3</td>
    <td>3</td>
    <td>3</td>
    <td>1</td>
    <td>3</td>
    <td>1</td>
  </tr>
  <tr>
    <th>Repeating Block</th>
    <td colspan="2">1</td>
    <td colspan="3">2</td>
    <td colspan="1">3</td>
    <td colspan="2">1</td>
    <td colspan="3">2</td>
    <td colspan="1">3</td>
    <td>1</td>
    <td colspan="3">2</td>
    <td colspan="1">3</td>
    <td colspan="1">2</td>
    <td colspan="1">3</td>
  </tr>
  <tr>
    <th>Repeating Unit</th>
    <td colspan="6">1</td>
    <td colspan="6">2</td>
    <td colspan="5">3</td>
    <td colspan="2">4</td>
  </tr>
</table>

### Random Copolymers
Random copolymers can be generated. The random sequence is generated by first creating a list where the particle type of each segment is given as for a block-like copolymer. This list is then shuffled using a [Knuth Shuffle](https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle). Now the particles are sequentially assigned to segments of its type. This way the particle types are randomly distributed, but the particles of one type are of increasing index along the chain. Note that if you want the exact sequence printed in your output file you need to set `itestpart = 10` in the `nmlParticle`.
* in `nmlParticle` set the parameters of the chain type `ict` which is to be random:

```fortran
txcopolymer(ict) = "random",
```

## Testin
The Testin-file covering these functionalities are `chain.sequence.fnw.mc.in` and `copoly-seq.mc.in`.

# Celllist {#Celllist}
## Overview
The cell list is implemented in molsim acts a a neighbor list, so that fewer particle pairs with a distance larger than \f$r_\mathrm{cut}\f$ are considered. This allows for much quicker simulations, especially for large systems with short ranged interactions. The simulation duration scales with \f$\mathcal{O}(n)\f$.
The advantages are compared to the other neighbor lists are the following:
 * The cell list is updated after every accepted step, so that is is always up to date. This allows for non-local moves to be carried out without needing to artificially enlarge the cell size. This is not the case for `vlist` and `llist` which are currently implemented in molsim.
 * The cell structure is generated only once at the beginning of the simulation. The time for the generation of the cell list scales with the number of cells, but the execution speed should be independent of the number of cells. In `llist` the number of cells were capped to a small number (1000), so that large systems had to use large cells. Increasing the number of possible cells did not help either, as the time to find the neighbors in each step increased with the number of cells, making the `llist` very slow when having many cells (\f$\approx 10^6\f$).


## Details of Implementation
### Working Procedure
The usage of a neighbor list within Molsim is controlled by the namelist `nmlIntList`. Usage if the cell list is achieved by setting ``txintlist = 'cellList'``. The cell size is the chosen to be as small as possible to fill the whole space evenly with cells, while still remaining larger than `rcut`.

One can create cells of different size, by changing the value of `drnlist`, which will be added to `rcut` for setting up the cell list. Most notably one can set `drnlist` to a negative value to even further reduce the cell size. This will lead to less particles with a distance larger than `rcut`, but more cells will be considered. Whether or not this approach increases the speed of the simulation may depend on the system.

All parameters except `txintlist` and `drnlist` of `nmlIntList` are not used for the cell list. A typical input file could be the following:

```fortran
&nmlIntList txintlist = 'cellList', drnlist = 0.0, /
```

## Testin
An example for the usage of the cell list can be found the `Testin` directory.
* `lj.mc.clist.in`

## Further Reading
* cell List https://en.wikipedia.org/wiki/Cell_lists

# SSO {#SSO}
## Overview
###   SSO
In this section the sso-functionality is added to Molsim. For a detailed description of the algorith see the master-thesis of pascal or the corresponding [publication](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00797). The SSO algorithm allows for on the fly optimization of the displacement parameter of the single particle move during the equilibration run.
### Sperarate Local from Non-Local
In addition the `lmcsep` parameter was introduced. This separates local (Single-Particle and SSO Move) from non-local moves (all the others) during the simulation. During each MC-Pass it is evaluated, wheter non-local moves are to be carried out (accoring to their probabiliy) and only if they will be carried out the neihbour list ist updated with a large drnlist (4* contou length of the longest chain). This allows for using a tigh neighbour-list during the simulation with local-moves.

## Testin
The following test-projects are added to the Testin directory:
* `sso.sso.in`: Test of the sso-move
* `sso.lmcsep.in`: Test of the lmcsep-setting

# Image generation with VMD {#VMD}

## Overview
In order to generate images for publications, presentations and for visual insight in simulated systems, Molsim offers the possibility to generate image files in the [VTF format ("VTF Trajectory Format")](https://www.github.com/olenz/vtfplugin/wiki/VTF-format). The VTF format has specifically been designed for people working with their own simulation code with the aim of visualization via the image software [VMD ("Visual Molecular Dynamics")](http://www.ks.uiuc.edu/Research/vmd). VMD offers a broad spectrum of functionality with regards to the visualization of molecular systems with the possibility of rendering high quality images.

## Details of Implementation
### File Types
The visualization of simulation snapshots requires two different file types.

* `vtf`: The vtf file contains atom, bond and box length definitions, followed by an arbitrary number of coordinate blocks.
* `tcl`: The tcl file forms the interface between Molsim and VMD. The whole settings of VMD may be controlled via tcl scripts. Here, color specifications, additional geometric objects (_e. g._ simulation cell) and the loading of the vtf file is stored.

### Working Procedure
The generation of images with Molsim is controlled by the namelists `nmlSystem`, `nmlImage` and namelist `nmlVTF`/`nmlVRML`. Note, that this manual covers only the image generation using the `vtf` format for the visualization with VMD. For the generation of images in the `vrml` format, please confer the MOLSIM manual.
In the following the work flow for the image generation is described.
* Request the generation of images by setting `limage = .true.` in the namelist `nmlSystem`

```fortran
&nmlSystem limage = .true.,/
```
* Request the generation of images using the `vtf` format by setting `lvtf = .true.` in the namelist `nmlImage`

```fortran
&nmlImage lvtf = .true.,/
```
* Further specifications may be given in the namelist `nmlVTF`. With just the default settings, images of the initial and of the final configuration will be stored, respectively.
* Run the simulation
* Open the `vtf` image file by executing the accompanying `tcl` script in VMD

```shell
vmd -e {name-of-project}.tcl
```
* Apply individual modifications to representation

### Render a snapshot
Here are some tipps on how to generate nice images using `VMD`.
* use the `AOShiny` material
* turn on the lights, often the scene generated looks less bright when rendered as presented one the screen in `VMD`

```
light 2 on
light 3 on
```
* turn on the shadows with `display shadows on`
* turn on the ambientocclusion `display ambientocclusion on`
* use the GLSL rendermode `display rendermode GLSL`

Then render the scene (with the path to you Tachyon executable adapted
```shell
render Tachyon scene.dat
exec "<path to you Tachyon executable>" -aasamples 48 scene.dat -res 2096 1048 -format BMP -o scene.bmp
```
You might want to adapt the name of the image (`scene.*`), the resolution and the the aasamples.

## Testin
An example for the usage of the image generation in the `vtf` format may be found in the `Testin` directory.
* vtfimage.fnw.mc.in

## Further Reading
* VTF [https://www.github.com/olenz/vtfplugin/wiki/VTF-format](https://www.github.com/olenz/vtfplugin/wiki/VTF-format)
* TCL [https://www.tcl.tk](https://www.tcl.tk)
* VMD [http://www.ks.uiuc.edu/Research/vmd](http://www.ks.uiuc.edu/Research/vmd)
