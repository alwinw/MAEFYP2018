# Bash cheatsheet 
  A; B    Run A and then B, regardless of success of A
  A && B  Run B if A succeeded
  A || B  Run B if A failed
  A &     Run A in background.

  ~% FILE="example.tar.gz"
  ~% echo "${FILE%%.*}"
  example
  ~% echo "${FILE%.*}"
  example.tar
  ~% echo "${FILE#*.}"
  tar.gz
  ~% echo "${FILE##*.}"
  gz
