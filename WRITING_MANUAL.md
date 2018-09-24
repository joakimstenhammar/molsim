Manual guide
============

This manual guide describes the generation of the manual MOLSIM and the working procedure of extending the existing manual.

## 1. Generating the manual with Doxygen
For the generation of the manual the open-source software [Doxygen](http://www.stack.nl/~dimitri/doxygen/) is used. It will generate the manual from special comments in the source code. Several outputs are possible. Current, for MOLSIM is only the Latex-Output supported, from which a `.pdf`-file is generated.

If you want to generate the manual of MOLSIM, go in the **doc**-directory with
```
cd doc
```
Then start the generation with
```
make doc
```
It will be generated a sub-directory named `latex` and the manual named `documentation.pdf` as a `.pdf`-file in the superior directory. In the `latex`-directory it is possible to edit  the file `doxygen.sty` for changes in the output formatting. If this sub-directory is not needed, it can be deleted with
```
make clean
```

## 2. For Contributors: How to document your new code for the manual
Doygen run through the code and recognize the comments, which the contributors has written for the manual. Therefor the comments has to be written in the right format. This is explained in section 2.2. Furthermore the manual can be extended by additional text, which is written in an external file. How to write text in an external file and to include in in the manual, is explained in the following subsection.

### 2.1 External files
External files are written in [Markdown](https://daringfireball.net/projects/markdown/). The files have got the endings `.md`. Markdown texts are converted into HTML, but in contrast to HTML, the plain text is good readable. In the following the Markdown syntax is shortly described.
#### 2.1.1 Markdown
##### headers
To write headers, use `#`'s. The number of `#`'s determines the level. Up to 6 levels are supported. An example:
```
# This is a level 1 header

### This is level 3 header
````
##### bold text
To write bold text, use two stars at the beginning and the end of the text. An example:
```
**bold text**
```
##### italic text
To write italic text, use one stars at the beginning and the end of the text. An example:
```
*bold text*
```
##### source code
For including source code the text should be wrapped in backticks(\`). For whole code blocks three backticks at the beginning and the end of the block should be used. If `fortran` is written behind the first three backticks, the code is highlighted in fortran-style. An example:
<pre><code>
```fortran
code
```
</code></pre>
For a code within one line only one backtick at the beginning and the end should be used.

##### links
or an inline link the link text is followed by a URL and an optional link title which together are enclosed in a set of regular parenthesis. The link title itself is surrounded by quotes.
Examples:
```
[The link text](http://example.net/)
[The link text](http://example.net/ "Link title")
[The link text](/relative/path/to/index.html "Link title")
[The link text](somefile.html)
```
##### images
The syntax for images is similar to that for links. The only difference is an additional ! before the link text. It would be the best, to save the image, which should be linked, in the directory `doc`.

##### tables
A table consists of a header line, a separator line, and at least one row line. Table columns are separated by the pipe (|) character. An example:
```
First Header  | Second Header
------------- | -------------
Content Cell  | Content Cell
Content Cell  | Content Cell
```
Only simple tables are supported. For complex tables HTML should be used.

#### 2.1.2 Include an external file to the manual
To include an external .md-file to the manual, the file has to be saved in the directory `doc`. Then the file `doxy` in the same directory has to be open. This file contains several flags, which determine the output. The flag `INPUT` is used to specify the input files. With
```
INPUT          += your_file.md
```
doxygen will consider the file. Doxygen handles the files in the order, in which they are in the file `doxy`.

### 2.2 Comment the code
This subsection describes how to comment the source code. Variables, namelists and subroutines are commented. Each has to be comment in a special way. But for all of them principally things are valid, which are described in the following subsections.

#### 2.2.1 General description of commenting
Each comment has to be started at the beginning of a line with
```
!>
```
For multi-line comments further lines has to start with
```
!!
```

#### 2.2.2 Special commands
Doxygen works with special commands, which effects several things, for example the reference to a variable. All commands starts with a ``\``. In the following the most important commands are listed.

##### \page
The descriptions of all variables, namelists, and subroutines start with the command `\page`. The sense of this is to have all variables in an own page-object, on which can be referenced. Furthermore several pages can be grouped in one big page.
After the \page-command the internal name of the thing, which should be documented, and then the real name are following. Usually both are the same and it is sufficient to write the real name once. But in cases, where the name of a variable appears in two or more namelists, the internal name should be write in this way:
```
namelistname_variablename
```
An example:
```
!> \page nmlSPDF_vtype vtype
```
The first line of all comments has this format.

##### \ref
With the command `\ref` it can be referenced to a page-object. If the name of a page-object appears in a description, it should be written the `\ref`-command ahead of the name.
An example:
```
!! The namelist  \ref nmlStatic contains variables that control the interval of the analysis and static analysis
```
##### \f$
Complex formulas in the .tex-format can be formulated with the `\f$`-command. The expression of the formula starts with `\f$` and ends with `\f$`
An example:
```
!! \f$ \langle r0\rangle \f$, and the projection of the molecular axes on the box axes denoted ``<x'>``, ``<y'>``, and ``<z'>``.
```
#### 2.2.3 Comment a variable
If the variable is declared in a module, the comments can be written directly above the declaration of the variable. If the declaration is within a subroutine, the comments has to be written above the subroutine.
The first line contains the `\page`-command. In the second line the type of the variable is enclosed in backticks. In the third line is the default value, if someone exists. At the beginning `default:` has to be enclosed in two stars and then the value in backticks. If the default value contains a variable, the variable should not be enclosed backticks, but a `\ref`-command should be written ahead of it. Two examples:
```
!! **default:** `0`
```
```
!! **default:** \ref nnwt*`0`
```
In the next line the variable is described. It starts with a star `*`. If there are multiple options for the value of the variable, all options should be described. The value has to be enclosed in backticks. A `*` has to be written ahead of all values. An example:
```
!> \page ldist
!! `logical`
!! **default:** `.false.`
!! * `.true.` Distribution functions are calculated (by routine DistFunc). Further specification is given in namelist \ref nmlDist.
!! * `.false.` No calculation of distribution functions.
```
#### 2.2.4 Comment a namelist
The comments of the namelists should be written at the beginning of the file, where the variables of the namelists are declared. The first line contains the `\page`-command. In the next lines a short description follows. Here the description should not start with a `*`. After the description the next line has to look so:
```
!! * Variables:
```
In the following lines the variables of the namelist are listed. For that the command `\subpage` is used. The command linked the variables to the namelist. Ahead of the `\subpage`-command a `*` has to be written, but the `*` has to be indent in comparison to the `*` in the line `!! * Variables:`. An example for a complete namelist:
```
!> \page nmlMSD
!! The namelist  \ref nmlMSD contains variables that control the calculation of the mean square displacement.
!! * Variables:
!!  * \subpage nmlMSD_sf
!!  * \subpage nmlMSD_cfin
```
To include the namelist in the manual, the namelist has to be include in the table of the file `input.md` in the directory `doc`.

#### 2.2.5 Comment a subroutine or module
Only subroutines, which are not within an other subroutine, should be commented. The first line contains the `\page`-command. The internal name of the page has to be the name of the file, which contains the subroutine or module, but without the ending `.F90`. The title has to be the complete filename. The second line contains the name of the subroutine/module enclosed in `**` and the third line contains the description enclosed in `*`. To distinguish the comments of subroutines/modules and variables easily, the comments of subroutines/modules should be enclosed of a common comment line with multiple stars. An example:
```
!************************************************************************
!> \page dynamic dynamic.F90
!! **DynamicDriver**
!! *driver of dynamic analysis routines*
!************************************************************************
```
#### 2.2.6 Comment a function
Functions are not commented, because they are within subroutines.
