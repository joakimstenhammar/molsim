Contributing guide
==================
This contributing guide describes the working procedure of the further development of Molsim. Changes and new features should always be implemented following the guidelines given below. The procedure can be summed up to the following steps:

  1. [Generate a new branch](#1-how-to-generate-a-new-branch)
  2. [Commit features/changes](#2-how-to-commit-featureschanges)
  3. [Coding Guidelines](#3-coding-guidelines)
  4. [Request a merge into Master](#4-how-to-request-a-merge-into-master)
  5. [Review a merge request](#5-how-to-review-a-merge-request)
  6. [Resolve all discussions](#6-how-to-resolve-all-discussions)
  7. [Merge it!](#7-how-to-merge-it)
  8. [Finalize your merge](#8-how-to-finalize-your-merge)

Most of the tasks may either be done using the [github interface](https://github.com/joakimstenhammar/molsim) or by using the command line. In this contributing guide the respective way is described, of which the authors think it is most convenient.

## 1. How to generate a new branch
Changes of Molsim may only be made within the scope of [branches](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging). This way you can implement your changes without messing with the main version (`master`) of Molsim. To generate a new branch:
```sh
git checkout -b <NAME-OF-NEW-BRANCH>
```
The just created branch exists only locally, to set it up on the remote repository use
```sh
git push --set-upstream origin <NAME-OF-NEW-BRANCH>
```
You're now in your new branch. You may now begin to implement new features/changes. Molsim comes with a feature called `Testin`. It allows for an efficient way to check whether changes of the Molsim source code lead to any malfunction of the program. This testing is an integral part of the extension of the Molsim functionality and needs to be performed at least once, before a branch may be merged into the master branch. Please follow the steps described in the the section *Testin* in the manual in order to run the `Testin` directory.

## 2. How to commit features/changes
The changes you apply to the code should be as efficient and non-invasive as possible. Try to divide your modifications in logically-associated chunks of code. These chunks can then individually be committed and described in a commit message. After you changed something you first have to stage the files in which changes were made and which you'd like to commit.
```sh
git add <LIST-OF-MODIFIED-FILES>
```
After staging one or several files you may now commit using
```sh
git commit -m "<YOUR-MESSAGE>"
```
The advantage of chopping all changes into smaller chunks of code is the option to revert individual commits. In addition, it enhances the transparency of what you do and others can better follow your changes. Whenever you stop working on your branch, you should always push your commits in order to upload them.
```sh
git push
```
Whenever you start to work on your branch again you should pull the branch in order to download the current version of your branch.
```sh
git pull
```
When you're done with your modifications, you may request to merge your branch into the `master`. Before doing so, please confer [this checklist](#appendix-checklist) and make sure, that all requirements have been satisfied.

## 3. Coding Guidelines

The Molsim code is not very well standardized and does not adhere to any standard. Just some best practices:
* When adding new subroutines, try keep them similar to existing subroutines.
* Indent blocks by 3 characters.
* Use genuine spaces rather than tabs.
* Try to confine your line width to 80 characters.
* Avoid putting multiple statements on the same line. It is not good for readability..
* If you add substantial parts, consider creating a new module for your changes. Place this module in a separate file.
* Use `END SUBROUTINE <name>` or `END FUNCTION <name>` instead of just `END`.
* Use IMPLICIT NONE in all program units.
* Avoid the GOTO statement.
* Avoid the use of ‘magic numbers’. A ‘magic number’ is a numeric constant that is hard wired into the code.
* Insert meaningful comments frequently.
* **Describe the Changes in the Manual**: A documentation on how to write the manual is given [here](WRITING_MANUAL.md).

## 4. How to request to merge into Master
In order to request a merge of your branch into the `master` browse the [gitub interface](https://github.com/joakimstenhammar/molsim/pulls) to create a new merge request. Note to request the merge with a WIP-prefix. The merge is then designated as "work in progress".

## 5. How to review a merge request
When you are assigned to a merge request you are supposed to read the changes in the code and check for possible mistakes. Besides of possible mistakes the code should be straight forward to understand. If parts of the code are difficult to understand, request more comments! When reviewing the code differences using git (or the [github interface](https://github.com/joakimstenhammar/molsim/pulls)) hide whitespace changes. Whenever you find something to be wrong or incomprehensible, add a comment to the respective line of code. This will start a discussion to be resolved by the author of the modification.

Some general rules when commenting the changes of the code:
* When reviewing the code, mark comments regarding ''cosmetic'' changes in the code with :sparkles: (`:sparkles:`).

Besides of the reviewing of the detailed code modifications, please check whether all requirements in [this checklist](#appendix-checklist) have been met. When you are finished reviewing all modifications and all discussions have been resolved, remove the WIP-prefix from the merge request. This signals the author of the modifications, that his branch may be merged.

## 6. How to resolve all discussions
When the assignee of your merge request has fully reviewed your modifications, it is your responsibility to resolve all discussions by marking them as resolved.

Discussions with a :sparkles: in it are only of cosmetic nature. You are not obliged to resolve this kind of discussions. Still, it is better to resolve it or discuss with the assignee, why for example it should not be changed.

## 7. How to merge it
When finally all discussions have been resolved and the WIP-prefix has been removed  by the assignee of the merge request, you may merge your branch into the `master`. But first you need to perform some housekeeping:

### 7.1 Update the Changelog
Update the [changelog](CHANGELOG.md) according to [this description](http://keepachangelog.com/en/0.3.0/). If you created backwards incompatible changes in the input, please describe your changes in the changelog.

### 7.2 Declare Stable
Declare your version as the current stable version. Do this by running `make declarestable` in the `Testin` directory and commit/push the new `stable.md5sum` file.

### 7.3 Merge Branch
Merge your branch by using the [github interface](https://github.com/joakimstenhammar/molsim/pulls).

## 8. How to finalize your merge
After the merge has been done, there are still a few things in order to finalize your merge request:

### 8.1 Generate a tag
Generating a tag means to mark this specific point of the version history. First checkout the `master` branch, and pull it:
```sh
git checkout master
git pull
```
Then tag the current commit:
```sh
git tag vX.Y.Z
```
where `X`, `Y` and `Z` correspond to the major (X), minor (Y) and patch level (Z) number of the software version ([Semantic Versioning](http://semver.org/)). 

After tagging
```sh
git push --tags
```
in order to synchronize your new tag.

## Appendix: Checklist
* [ ] Make sure, that the `Testin` runs clean! For further information confer the [manual](documentation.pdf).
* [ ] Ascertain, that the compilation with `mode=warn` does not trigger any warnings
* [ ] ascertain, that running a test input file covering your modifications runs clean after compiling with `mode=debug`
* [ ] Make sure that your code runs in parallel mode, and if not add an error message which prevents from running your code in parallel.
* [ ] Describe all changes in the manual ([how to write the manual](WRITING_MANUAL.md))
* [ ] Change the version number corresponding to [Semantic Versioning](http://semver.org/)!
