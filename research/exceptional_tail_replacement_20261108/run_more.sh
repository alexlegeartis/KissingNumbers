set -e
for d in 27 29 30 31; do
 OPENBLAS_NUM_THREADS=3 OMP_NUM_THREADS=3 taskset -c 0-3 python research/exceptional_tail_replacement_20261108/global_tetrads.py "$d" --seconds 90 --out "research/exceptional_tail_replacement_20261108/tetrads$d" > "research/exceptional_tail_replacement_20261108/tetrads$d.log" 2>&1
done
