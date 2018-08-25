git filter-branch --force --index-filter \
  "git rm -rf --cached --ignore-unmatch session-files" \
  --prune-empty -- --all
# git for-each-ref --format="%(refname)" refs/original/ | xargs -n 1 git update-ref -d
git for-each-ref --format='delete %(refname)' refs/original | git update-ref --stdin
git gc

git filter-branch --force --index-filter \
  "git rm -rf --cached --ignore-unmatch R/.RData" \
  --prune-empty -- --all
# git for-each-ref --format="%(refname)" refs/original/ | xargs -n 1 git update-ref -d
git for-each-ref --format='delete %(refname)' refs/original | git update-ref --stdin
git gc

git filter-branch --force --index-filter \
  "git rm -rf --cached --ignore-unmatch *.lyx#" \
  --prune-empty -- --all
# git for-each-ref --format="%(refname)" refs/original/ | xargs -n 1 git update-ref -d
git for-each-ref --format='delete %(refname)' refs/original | git update-ref --stdin
git gc

git reflog expire --expire=now --all
git gc --aggressive --prune=now
