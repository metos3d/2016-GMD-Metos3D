#!/bin/env python

header = '''#!/bin/bash
#PBS -T intmpi
#PBS -b 8
#PBS -l cpunum_job=24
#PBS -l elapstim_req=24:00:00
#PBS -l memsz_job=10gb
#PBS -N speedtest
#PBS -o speedtest.out
#PBS -e speedtest.out
#PBS -q clfo2

cd $PBS_O_WORKDIR

. ../../petsc/petsc-3.4.5.env
'''

mpicmd = '''mpirun $NQSII_MPIOPTS -n %d ./tmmmitgchemdic -numtracers 6 -i dicini.petsc,alkini.petsc,po4ini.petsc,dopini.petsc,o2ini.petsc,feini.petsc \
	-me Ae -mi Ai -t0 0.0 -iter0 0 -deltat_clock 0.0013888888888889 -max_steps 720 -write_steps 7200000 -o dic.petsc,alk.petsc,po4.petsc,dop.petsc,o2.petsc,fe.petsc \
	-external_forcing -use_profiles -biogeochem_deltat 43200.0 -periodic_matrix -matrix_cycle_period 1.0 -matrix_cycle_step 0.0833333333333333 \
	-periodic_biogeochem_forcing -periodic_biogeochem_cycle_period 1.0 -periodic_biogeochem_cycle_step 0.0833333333333333 '''

print header
for i in range(0,192):
	for j in (1,2,3,4,5):
		print mpicmd % (i+1)

