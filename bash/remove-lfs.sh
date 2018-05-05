/bin/bash

#============================
# Remove LFS from Repo
# Alwin Wang
#----------------------------
# https://stackoverflow.com/questions/48699293/how-to-i-disable-git-lfs
# https://github.com/git-lfs/git-lfs/issues/910#issuecomment-238389388

# git lfs ls-files | cut -d ' ' -f 3 > lfs-files.txt
# cat bash/lfs-files.txt | xargs touch

# git rm -r --cached .
# git add .

git lfs uninstall
git filter-branch -f --prune-empty --tree-filter '
  git lfs checkout
  git lfs ls-files | cut -d " " -f 3 | xargs touch
  git rm -f .gitattributes
  git lfs ls-files | cut -d " " -f 3 | git add
' --tag-name-filter cat -- --all


git filter-branch --prune-empty --tree-filter '
git config -f .gitconfig lfs.url "http://127.0.0.1:8080/user/repo"
git lfs untrack "*.npy"
git add .gitattributes .gitconfig

for file in $(git ls-files | xargs git check-attr filter | grep "filter: lfs" | sed -r "s/(.*): filter: lfs/\1/"); do
  echo "Processing ${file}"

  git rm -f --cached ${file}
  echo "Unadding $file lfs style"
  git add ${file}
done' --tag-name-filter cat -- --all

git filter-branch --prune-empty --index-filter 'git rm --ignore-unmatch --cached "*.mov"'


git filter-branch --prune-empty --tree-filter '
git lfs untrack "*.dump"
git add .gitattributes .gitconfig

for file in $(git ls-files | xargs git check-attr filter | grep "filter: lfs" | sed -r "s/(.*): filter: lfs/\1/"); do
  echo "Processing ${file}"

  git rm -f --cached ${file}
  echo "Unadding $file lfs style"
  git add ${file}
done' --tag-name-filter cat -- --all

