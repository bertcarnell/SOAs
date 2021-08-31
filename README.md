# SOAs

Creates Stratum Orthogonal Arrays (a.k.a. Strong OAs)

- **Author**:  Ulrike Groemping
- **Contributor**:  Rob Carnell

|Actions|
|:-----:|
|[![R-CMD-check](https://github.com/bertcarnell/SOAs/actions/workflows/r_cmd_check.yml/badge.svg)](https://github.com/bertcarnell/SOAs/actions/workflows/r_cmd_check.yml)|

## System Dependencies

For the `arrangements` package:

- Ubuntu: `apt-get update; apt-get install libgmp3-dev -y`
- MacOS: `brew install gmp`

## Development Process

RStudio is highly recommended to make package development easier.

1. Open the RStudio project `SOAs.Rproj`
2. Pull from github
    - Resolve any issues
3. Make changes
    - If the package functions are needed during development, use Load All
        - Build tab -> More -> Load All (ctrl-shift-l)
        - OR `devtools::load_all()`
    - Writing new tests
        - Open the `.R` file in RStudio that you would like to create a test for
        - `usethis::use_test()`
    - Inserting Roxygen comments
        - In RStudio, position your cursor in the function to be documented
        - Code -> Insert Roxygen Skeleton
4. **Document** using `roxygen2`:
    - Build tab -> More -> Document (ctrl-shift-d)
    - OR `devtools::document()`
5. **Test** using `testthat`:
    - Build tab -> More -> Test (ctrl-shift-t)
    - OR `devtools::test()`
6. **Check**:
    - Build tab -> Check
    - OR `rcmdcheck::rcmdcheck()`
    - OR `R CMD build SOAs; R CMD check SOAs_x.y.z.tar.gz`
7. (If desired) **Install** the package locally
    - Build tab -> Install and Restart
8. **Git** Process
    - Pull - to get the latest - `git pull`
    - Add - changed files to be committed - `git add`
    - Commit - to the local repo - `git commit -m "message"`
    - Push - to the remote repo - `git push`
