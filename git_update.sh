git checkout dev/bash_runs          &&
git merge master                    &&
git checkout dev/cpp_semtex         &&
git merge master                    &&
git checkout dev/R_airfoil-analysis &&
git merge master                    &&
git checkout dev/R_batch-analysis   &&
git merge master                    &&
git checkout dev/R_single-analysis  &&
git merge master                    &&
git checkout wip/doc_reports        &&
git merge master                    &&
git checkout master                 &&
git push --all
