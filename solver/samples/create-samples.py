#!/usr/bin/env python

import os
import sys
import numpy as np
import numpy.random as npr

#
#   read_PETSc_vec
#
def read_PETSc_vec(filepath):
	# debug
	print "reading in PETSc vector ...", filepath
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
#   write_PETSc_vec
#
def write_PETSc_vec(v, filepath):
	# debug
	print "writing PETSc vector ...", filepath
	# open file
	f = open(filepath, "wb")
	# write PETSc vector header
	np.array(1211214, dtype = ">i4").tofile(f)
	# write length
	np.array(len(v), dtype = ">i4").tofile(f)
	# write values
	np.array(v, dtype = ">f8").tofile(f)
	# close file
	f.close()

#
#   read_PETSc_mat
#
def read_PETSc_mat(filepath):
	# debug
	print "reading in PETSc matrix ...", filepath
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
#	create_samples_read_data
#
def create_samples_read_data():
	# debug
	print "read data ..."
	# set file paths
	lsmpath = m3dprefix + "/data/data/TMM/2.8/Geometry/landSeaMask.petsc"
	volpath = m3dprefix + "/data/data/TMM/2.8/Geometry/volumes.petsc"
	# read in land sea mask and volumes
	# convert lsm to a matrix of ints
	lsm = read_PETSc_mat(lsmpath)
	lsm = lsm.astype(int)	
	vol = read_PETSc_vec(volpath)
	# return data
	return lsm, vol

#
#	create_samples
#
def create_samples(m3dprefix, argv):
	# debug
	print "create samples ..."
	# retrieve mass
	mass = float(sys.argv[1])
	# debug
	print "using mass ...", mass, "mmol P m^(-3)"
	# retrieve tracers
	tracers = sys.argv[2:]
	ntracer = len(tracers)	
	# debug
	print "using tracer(s) ...", tracers
	# read data
	lsm, vol = create_samples_read_data()
	# normalize volumes
	volnorm = vol / sum(vol)
	# divide mass by no of tracers
	masstracer = mass / ntracer
	# loop over tracers
	for itracer in range(ntracer):
		# get random values
		val = 0.5 + 0.25 * npr.ranf(vol.shape)
		# scale
		val = masstracer * val / sum(val * volnorm)
		# debug
# 		print val
# 		print min(val), max(val)
# 		print sum(val * volnorm)
		# save vector
		filepath = "init.%s.petsc" % tracers[itracer]
		write_PETSc_vec(val, filepath)

#
#   main
#
if __name__ == "__main__":
	# no arguments?
	if len(sys.argv) <= 2:
		# print usage and exit with code 1
		print("usage: ./create-samples.py [MASS...] [TRACERS...]")
		print("example:")
		print("  $> ./create-samples.py 2.17 N P Z D DOP")
		sys.exit(1)
	# determine path prefix for geometry data
	m3dprefix = os.path.expanduser("~/.metos3d")
	# create samples
	create_samples(m3dprefix, sys.argv)
