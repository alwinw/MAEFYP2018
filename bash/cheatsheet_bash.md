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


Sample task
task(){
   sleep 0.5; echo "$1";
}
Sequential runs
for thing in a b c d e f g; do 
   task "$thing"
done
Parallel runs
for thing in a b c d e f g; do 
  task "$thing" &
done
Parallel runs in N-process batches
N=4
(
for thing in a b c d e f g; do 
   ((i=i%N)); ((i++==0)) && wait
   task "$thing" & 
done
)
It's also possible to use FIFOs as semaphores and use them to ensure that new processes are spawned as soon as possible and that no more than N processes runs at the same time. But it requires more code.

N processes with a FIFO-based semaphore:
open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}
run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
    "$@" 
    printf '%.3d' $? >&3
    )&
}

N=4
open_sem $N
for thing in {a..g}; do
    run_with_lock task $thing
done 
