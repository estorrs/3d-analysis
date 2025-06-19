export PATH="/miniconda/bin:$PATH"
export LSF_DOCKER_VOLUMES="/storage1/fs1/dinglab:/storage1/fs1/dinglab /scratch1/fs1/dinglab:/scratch1/fs1/dinglab /home/estorrs:/home/estorrs"
LSF_DOCKER_PORTS='8181:8888' bsub -R 'select[port8181=1] rusage[mem=100GB] span[hosts=1]' -M 100GB -q general -G compute-dinglab -oo log.txt -a 'docker(estorrs/3d_analysis:latest)' 'jupyter notebook --port 8888 --no-browser --ip=0.0.0.0'

# log job output
# bsub -G ${group_name} -I -q general -a 'docker_exec(<JOBID>)' cat /path/to/log/file
