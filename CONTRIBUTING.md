# Contribution Guidelines

## Reporting issues

- **Search for existing issues.** Please check to see if someone else has reported the same issue.
- **Share as much information as possible.** Include operating system and version, browser and version. Also, include steps to reproduce the bug.

## Project Setup
Refer to the [README](README.md).

## Code Style

### Variable Naming
Not all current code follows the conventions below but these will be followed for future developments. 
- Maximize the use  of semantic and descriptive variables names (e.g. `faceIndices` not `fcInd` or `fi`). Avoid abbreviations except in cases of industry wide usage. In some cases non-descriptive and short variable names are exceptable for instance vertices (points), faces, edges, colors and logic arrays may be denoted `V`, `F`, `E`, `C`, `L`. Furthermore, if a mathematrical symbol or letter is commonly used for some entity it may be acceptable to use short names e.g. coordinates may be referred to as `X`, `Y` and `Z` and image coordinates of indices may be referred to as `I`, `J` and `K`. 

## Testing
Proving example files as well as a description of the expected outcomes. This can be in the form of a log file or data file that contain partial solution data. 

## Pull requests
- Try not to pollute your pull request with unintended changes – keep them simple and small. If possible, squash your commits.
- Try to share how your code has been tested before submitting a pull request.
- If your PR resolves an issue, include **closes #ISSUE_NUMBER** in your commit message (or a [synonym](https://help.github.com/articles/closing-issues-via-commit-messages)).
- Review
    - If your PR is ready for review, another contributor will be assigned to review your PR
    - The reviewer will accept or comment on the PR. 
    - If needed address the comments left by the reviewer. Once you're ready to continue the review, ping the reviewer in a comment.
    - Once accepted your code will be merged to `master`
