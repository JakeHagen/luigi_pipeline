
#!/bin/bash
#BSUB -J luigi_star_align_test 
#BSUB -q expressalloc 
#BSUB -P acc_argmac01a
#BSUB -n 6  
#BSUB -R rusage[mem=8000]
#BSUB -R span[hosts=1]
#BSUB -W 1:00
#BSUB -o %J.stdout
#BSUB -eo %J.stderr



luigid --background --pidfile /hpc/users/hagenj02/luigi-server/luigid.pid --logdir /hpc/users/hagenj02/luigi-server/logs/ --state-path /hpc/users/hagenj02/luigi-server/state.pickler
python /hpc/users/hagenj02/luigi_pipeline/rnaseq_analysis.py all_star_align --workers 6


