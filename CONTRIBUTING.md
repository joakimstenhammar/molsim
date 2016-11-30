Contributing guide
==================
This contributing guide describes the working procedure of the further development of Molsim. Changes and new features should always be implemented following the guidelines given below. The procedure can be summed up to the following steps: 
 
  1. Generate a new branch  
  2. Commit features/changes  
  3. Request a merge into Master  
  4. Review a merge request  
  5. Resolve all discussions  
  6. Merge it!  
  7. Finalize your merge  

Most of the tasks may either be done using a browser ([https://git.rwth-aachen.de](https://git.rwth-aachen.de)) or by using the command line. In this contributing guide the respective way is described, of which the authors think it is most convenient.

## How to generate a new branch
Changes of Molsim may only be made within the scope of an isolated branch. This way the main version (`Master`) of Molsim is less probably affected by potential flaws. In order to generate a new branch browse [https://git.rwth-aachen.de](https://git.rwth-aachen.de) and access the project Molsim. Here you can find the tab `Branches`, where in turn you find the button `New Branch`. Using this button will lead to you being asked for a `Branch Name`. Further you'll be asked which other branch shall serve as the source for your new branch. Usually the `Master Branch` should be used as the parent version of your new branch. Finally use the `Create branch` button. You can now access your new branch by entering your git-directory and
```sh
git pull
git checkout NAME-OF-YOUR-BRANCH
```
You're now in the directory of your new branch. You may now begin to implement new features/changes.
## How to commit features/changes
The changes you apply to the code should be as efficient and non-invasive as possible. Try to divide your modifications in logically-associated chunks of code. These chunks can then individually be commited and described in a commit message. After you changed something you first have to `stage` the files in which changes were made and which you'd like to commit.
```sh
git add LIST-OF-MODIFIED-FILES
``` 
After staging one or several files you may now commit using
```sh
git commit
```
A vim-instance will open. Here you should enter a one-line description of the changes coming with this commit. When you're done with the description, save and close the vim-instance (`:wq`/`:x`). The advantage of chopping all changes into smaller chunks of code is the option to revert individual commits. In addition it enhances the transparency of what you do and others can better follow your changes. Whenever you stop working on your branch, you should always push your commits in order to upload them.
```sh
git push
```
Whenever you start to work on your branch again you should pull the branch in order to download the current version of your branch.
```sh
git pull
```
When you're done with your modifications, you may request a merge of your branch into the `master`.
## How to request a merge into Master
In order to request a merge of your branch into the `master` browse [https://git.rwth-aachen.de](https://git.rwth-aachen.de). Here you will find the tab `Merge Request`, where in turn you can use the `New Merge Request` button.
## How to review a merge request
When you are assigned to a merge request you are supposed to read the changes in the code and check for possible mistakes. The following tools can help you detect errors:
* run the tests in the `Testin` directory. The diff should yield no changes. If any file in the `Testin/save` or `Testin/in` dir was changed, review whether these are expected changes
* compile with `mode=warn`. Check that no warnings are affecting the parts of the code relevant for the merge request
* when viewing the code differences using git (or the gitlab interface) hide whitespace changes 

Some general rules when commenting the changes of the code:
* When reviewing the code mark comments regarding ''cosmetic'' changes in the code with :sparkles: (`:sparkles: `).
* Mark all general comments which require further attention with :negative_squared_cross_mark: (`:negative_squared_cross_mark:`). When the issue of the comment was resolved, mark it with :white_check_mark: (`:white_check_mark:`)

Some other things to check:
* are all new features described in the wiki?
* is the version changed (correctly? adhere to [Semantic Versioning](http://semver.org/))
## How to resolve all discussions
## How to merge it
## How to finalize your merge
