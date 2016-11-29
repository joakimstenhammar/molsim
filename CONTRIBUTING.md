Contributing guide
==================

## How to review a merge request

When you are assigned to a merge request you are supposed to read the changes in the code and check for possible mistakes. The following tools can help you detect errors:
* run the tests in the `Testin` directory. The diff should yield no changes. If any file in the `Testin/save` or `Testin/in` dir was changed, review wheteher these are expected changes
* compile with `mode=warn`. Check that no warnings are affecting the parts of the code relevant for the merge request
* when viewing the code differences using git (or the gitlab interface) hide whitespace changes 

Some general rules when commenting the code:
* When reviewing the code mark comments regarding ''cosmetic'' changes in the code with :sparkles: (`:sparkles: `).
* Mark all general comments which requre further attention with :negative_squared_cross_mark: (`:negative_squared_cross_mark:`). When the isse of the comment was resolved, mark it with :white_check_mark: (`:white_check_mark:`)