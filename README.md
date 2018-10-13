## MAE FYP 2018
**Examination of the details of 2D Vorticity Generation Around the Airfoil During Starting and Stopping Phases**

**By: Alwin Wang**

**Supervisor: Hugh Blackburn**

**Git Branches**
* `master`              : Stable changes, Report writing
* `afmc`                : Legacy only
  + `R_airfoil-analysis`  : Determing spline length, remeshing the airfoil, etc 
  + `R_batch-analysis`    : Batch load and process sesh, mesh and dump
  + `R_single-analysis`   : Single analysis and testing
  + `bash_runs`           : Bash scripts to run analysis
  + `cpp_semtex`          : Semtex development work
  + `doc_reports`         : Report and paper writing and research
* `dev`                 : Development work for code
  + `example`             : For examples (tbc)
  + `mesh`                : Generating new meshes, refining mesh and reading mesh
  + `semtex`              : Semtex development work
  + `unit-tests`          : Various test cases to test code
* `wip`                 : Work in progress for objectives
  + `convg`               : Convergence analysis
  + `docs`                : AFMC, FYP and lit review
  + `kutta`               : Establishment of Kutta condition
  + `vortgen`             : Detailed analysis of the vorticity generation

**Folder Structure**
* `analysis`            : Analyses conducted
  + `convg`               : Convergence analysis (RProj, sesh, bash)
    - `results`             : Semtex output of session file
    - `output`              : Output of R scripts (RData, images)
  + `panel`               : Panel method analysis (RProj, sesh, bash)
    - `results`             : Semtex output of session file
    - `output`              : Output of R scripts (RData, images)
  + `kutta`               : Kutta condition analysis (RPoj, sesh, bash)
    - `results`             : Semtex output of session file
    - `output`              : Output of R scripts (RData, images)
  + `vortgen`             : Vorticity generation analysis (RPoj, sesh, bash)
    - `results`             : Semtex output of session file
    - `output`              : Output of R scripts (RData, images)
* `docs`                : Documents including papers  (max pages/words indicated)
  + `afmc_abstract`       : Abstract for 21st AFMC    ( 1 page )
  + `afmc_full`           : Full paper for 21st AFMC  ( 4 pages)
  + `fyp1_proposal`       : Proposal for FYP          ( 5 pages)
  + `fyp2_progress`       : Progress report for FYP   (10 pages)
  + `fyp3_final`          : Final FYP report          (50 pages/20,000 words)
  + `fyp4_paper`          : Research paper            ( 6 pages)
  + `litr_papers`         : Literature review, papers, summaries, etc
  + `litr_textbooks`      : Selected textbooks, summaries, etc
  + `semtex`              : Semtex user guides, etc
* `git`                 : Scripts and other git resources
* `R`                   : Generic R scripts for multiple applications (RProj)
* `semtex`              : Semtex development work
* `tests`               : Session file, bash and results folder
  + `kovas`               : Kovasznay flow examples inc. convergence
  + `sp5`                 : 0 Original sp5 files
  + `naca0012`            : 1 Joined bndry_prf 
  + `remesh`              : 2 LE and TE clustered mesh
  + `panel`               : Testing of panel method
  
   
