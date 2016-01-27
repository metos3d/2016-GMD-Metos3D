#### Profiling experiments

Setup PETSc environment:

	$> . petsc-3.3-p4.env.txt

Compile model and run experiment (submit job):
	
- ``MITgcm-PO4-DOP`` model:
	
		$> cd MITgcm-PO4-DOP/
		$> metos3d simpack MITgcm-PO4-DOP
		$> qsub job/profile.MITgcm-PO4-DOP.job
	
- ``N`` model:
	
		$> cd N/
		$> metos3d simpack N
		$> qsub job/profile.N.job
	
- ``N-DOP`` model:
	
		$> cd N-DOP/
		$> metos3d simpack N-DOP
		$> qsub job/profile.N-DOP.job
	
...
	
- ``NPZD-DOP`` model:
	
		$> cd NPZD-DOP/
		$> metos3d simpack NPZD-DOP
		$> qsub job/profile.NPZD-DOP.job

Create figures:

- ``MITgcm-PO4-DOP`` model:

		$> ./create-figures.py MITgcm-PO4-DOP/

- ``N`` model:

		$> ./create-figures.py N/

...

- ``NPZD-DOP`` model:

		$> ./create-figures.py NPZD-DOP/
