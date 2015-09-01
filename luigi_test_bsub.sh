
#!/bin/bash
#BSUB -J luigi_star_align_test 
#BSUB -q low 
#BSUB -P acc_argmac01a
#BSUB -n 12  
#BSUB -R rusage[mem=3500]
#BSUB -R span[hosts=1]
#BSUB -W 4:00
#BSUB -o %J.stdout
#BSUB -eo %J.stderr



luigid --background --pidfile /hpc/users/hagenj02/luigi-server/luigid.pid --logdir /hpc/users/hagenj02/luigi-server/logs/ --state-path /hpc/users/hagenj02/luigi-server/state.pickler
python /hpc/users/hagenj02/luigi_pipeline/rnaseq_analysis.py all_star_align

exit


