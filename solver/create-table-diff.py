#!/usr/bin/env python -W ignore

import os
import sys
import numpy as np

# set model list
modellist = ['MITgcm-PO4-DOP', 'N', 'N-DOP', 'NP-DOP', 'NPZ-DOP', 'NPZD-DOP']
tracerlist = [['PO4', 'DOP'], ['N'], ['N', 'DOP'], ['N', 'P', 'DOP'], ['N', 'P', 'Z', 'DOP'], ['N', 'P', 'Z', 'D', 'DOP']]	

#
#   read_PETSc_vec
#
def read_PETSc_vec(filepath):
	# debug
# 	print "reading in PETSc vector ...", filepath
	# open file
	f = open(filepath, "rb")
	# omit header
	np.fromfile(f, dtype = ">i4", count = 1)
	# read length
	nvec = np.fromfile(f, dtype = ">i4", count = 1)
	# read values
	v = np.fromfile(f, dtype = ">f8", count = nvec)
	# close file
	f.close()
	# return vector
	return v

#
#   read_PETSc_mat
#
def read_PETSc_mat(filepath):
	# debug
# 	print "reading in PETSc matrix ...", filepath
	# open file
	f = open(filepath, "rb")
	# omit header
	np.fromfile(f, dtype = ">i4", count = 1)
	# read dims
	nx     = np.fromfile(f, dtype = ">i4", count = 1)
	ny     = np.fromfile(f, dtype = ">i4", count = 1)
	# read no of non zeros
	nnz    = np.fromfile(f, dtype = ">i4", count = 1)
	# non zeros per row
	nrow   = np.fromfile(f, dtype = ">i4", count = nx)
	# column indices
	colidx = np.fromfile(f, dtype = ">i4", count = nnz)
	# values
	val    = np.fromfile(f, dtype = ">f8", count = nnz)
	# close file
	f.close()
	# create full matrix (default np.float64)
	mat = np.zeros(shape = (nx, ny))
	offset = 0
	for i in range(nx):
		if not nrow[i] == 0.0:
			for j in range(nrow[i]):
				mat[i, colidx[offset]] = val[offset]
				offset = offset + 1
	# return matrix
	return mat

#
#	create table
#
def create_table(vlist, vol):
	# debug
# 	print("create table ...")
	print(r'''
%tab:norms
\begin{table*}
\setlength{\tabcolsep}{10pt}
\caption{\label{tab:norms}
Difference in Euclidean norm ...}
\vspace{0.3cm}
\small
\centering
\begin{tabular}{l S S}
\tophline
Model			&	{$\| \mathbf{y}_{\mathrm{SP}} - \mathbf{y}_{\mathrm{NK}} \|_2$}	& 	{$\| \mathbf{y}_{\mathrm{SP}} - \mathbf{y}_{\mathrm{NK}} \|_{2, \Omega}$}	\\
\middlehline''')
	# go through models
	for i in range(0, len(modellist)):
		# combine volumes and tracer vectors
		vols = np.tile(vol, len(tracerlist[i]))
		tracer = np.concatenate(vlist[i])
		# debug
		print("%-15s" % modellist[i]),
		print(r'''& '''),
# 		print(r'''& \num{'''),
# 		print("vector lengths: %d, %d" % (len(vols), len(tracer)))
		# compute 2-norm
		norm = np.linalg.norm(tracer)
		# debug
		print("%.3e" % norm),
		print(r''' & '''),
# 		print(r'''} & \num{'''),
		# compute weighted norm
		normw = np.linalg.norm(tracer * np.sqrt(vols))
		# debug
		print("%.3e" % normw),
		print(r''' \\''')
# 		print(r'''} \\''')
	print(r'''\bottomhline
\end{tabular}
\end{table*}
''')

#
#   main
#
if __name__ == "__main__":	
	# debug
# 	print("comparing models ... " + str.join(', ', modellist))
	# set prefix for newton experiments
	if len(sys.argv) == len(modellist)+1:
		prefixlist = sys.argv[1:len(modellist)+1]
	else:
		prefixlist = ('newton. '*len(modellist)).split(' ')
		prefixlist.pop()
	# debug
# 	print("using given prefix for model ... %s" % prefixlist)

	# read in data
# 	print("read in data ...")
	# diffs
	vlist = []
	for i in range(0, len(modellist)):
		modeldir = modellist[i]
# 		print("model: %s tracer: %s" % (modeldir, str.join(', ', tracerlist[i])))
		# add new empty list
		vlist.append([])
		for j in range(0, len(tracerlist[i])):
			tracer = tracerlist[i][j]
			# read vectors
			vsp = read_PETSc_vec("%s/work/spinup.%s.petsc" % (modeldir, tracer))
			vnk = read_PETSc_vec("%s/work/%s%s.petsc" % (modeldir, prefixlist[i], tracer))
			# compute diff
			vdiff = vsp - vnk
			# add to list
			vlist[i].append(vdiff)
	# volumes
# 	print("read in volumes ...")
	# determine path prefix for geometry data
	m3dprefix = os.path.expanduser("~/.metos3d")
	volpath = m3dprefix + "/data/data/TMM/2.8/Geometry/volumes.petsc"
	vol = read_PETSc_vec(volpath)
	# create table
	create_table(vlist, vol)
	
