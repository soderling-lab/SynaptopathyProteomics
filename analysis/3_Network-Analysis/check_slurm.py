#!/usr/bin/env python3

import os
vars = ['SLURM_JOBID','SLURM_CPUS_PER_TASK']
envars = [os.environ.get(var) for var in vars]


