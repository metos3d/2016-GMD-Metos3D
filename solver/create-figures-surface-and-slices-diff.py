#!/usr/bin/env python -W ignore

import os
import sys
import re
import numpy as np

import matplotlib
matplotlib.rc("font", **{"family" : "sans-serif"})
matplotlib.rcParams.update({'font.size': 14})
matplotlib.rc("text", usetex = True)
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp
from mpl_toolkits.basemap import Basemap

#
#   read_PETSc_vec
#
def read_PETSc_vec(file):
    # open file
    # omit header
    # read length
    # read values
    # close file
    f = open(file, "rb")
    np.fromfile(f, dtype=">i4", count=1)
    nvec = np.fromfile(f, dtype=">i4", count=1)
    v = np.fromfile(f, dtype=">f8", count=nvec)
    f.close()

    return v

#
#   read_PETSc_mat
#
def read_PETSc_mat(file):
    # open file
    f = open(file, "rb")
    # omit header
    np.fromfile(f, dtype=">i4", count=1)
    # read dims
    nx     = np.fromfile(f, dtype=">i4", count=1)
    ny     = np.fromfile(f, dtype=">i4", count=1)
    nnz    = np.fromfile(f, dtype=">i4", count=1)
    nrow   = np.fromfile(f, dtype=">i4", count=nx)
    colidx = np.fromfile(f, dtype=">i4", count=nnz)
    val    = np.fromfile(f, dtype=">f8", count=nnz)
    
#    print "dims"
#    print nx, ny
#    print "nnz"
#    print nnz
#    print "nrow"
#    print nrow
#    print "colidx"
#    print colidx
#    print "val"
#    print val
    
    # close file
    f.close()
    
    # create full matrix
    matfull = np.zeros(shape=(nx, ny), dtype = float)
    offset = 0
    for i in range(nx):
        if not nrow[i] == 0.0:
            for j in range(nrow[i]):
                matfull[i, colidx[offset]] = val[offset]
                offset = offset + 1

    return matfull

#
#   read_data
#
def read_data(m3dprefix, modeldir, prefix, tracer):
    # debug
	print("reading data ... %s" % tracer)
    
    # v1d, z, dz, lsm (land sea mask)
	v1dsp = read_PETSc_vec("%s/work/spinup.%s.petsc" % (modeldir, tracer))
	v1dnk = read_PETSc_vec("%s/work/%s%s.petsc" % (modeldir, prefix, tracer))
    
	v1d = np.fabs(v1dsp - v1dnk)

	z   = read_PETSc_vec("%s/data/data/TMM/2.8/Forcing/DomainCondition/z.petsc" % m3dprefix)
	dz  = read_PETSc_vec("%s/data/data/TMM/2.8/Forcing/DomainCondition/dz.petsc" % m3dprefix)
	lsm = read_PETSc_mat("%s/data/data/TMM/2.8/Geometry/landSeaMask.petsc" % m3dprefix)
	lsm = lsm.astype(int)
    
	# dims
	nx, ny = lsm.shape
	nz = 15
	# v3d
	v3d = np.zeros(shape=(3, nx, ny, nz), dtype = float)
	v3d[:,:,:,:] = np.nan
    
	# v1d -> (v3d, z, dz)
	offset = 0
	for ix in range(nx):
		for iy in range(ny):
			length = lsm[ix, iy]
			if not length == 0:
				v3d[0, ix, iy, 0:length] = v1d[offset:offset+length]
				v3d[1, ix, iy, 0:length] = z[offset:offset+length]
				v3d[2, ix, iy, 0:length] = dz[offset:offset+length]
				offset = offset + length

	return v3d

#
#   create_figure_surface
#
def create_figure_surface(figid, aspect, xx, yy, cmin, cmax, levels, slices, v3d):
    print "# Creating surface ..."
    
    # data
    vv = v3d[0,:,:,0]
    # shift
    vv = np.roll(vv, 64, axis=1)

    # plot surface
    plt.figure(figid)
    # colormap
    cmap = plt.cm.bone_r
    
    # contour fill
    p1 = plt.contourf(xx, yy, vv, cmap=cmap, levels=levels, origin="lower")
    # contour lines
    p2 = plt.contour(p1, levels=levels, linewidths = (1,), colors="k")

    # slices
    s1 = xx[np.mod(slices[0]+64, 128)]
    s2 = xx[np.mod(slices[1]+64, 128)]
    s3 = xx[np.mod(slices[2]+64, 128)]
#    print s1, s2, s3
    plt.vlines([s1, s2, s3], -90, 90, color='k', linestyles='--')
    # set aspect ratio of axes
    plt.gca().set_aspect(aspect)

    # colorbar
    cbar = plt.colorbar(p1, format = '%.1e', fraction = 0.024, pad = 0.03)

    # basemap
    m = Basemap(projection="cyl")
    m.drawcoastlines(linewidth = 0.5)

    # xticks
    plt.xticks(range(-180, 181, 45), range(-180, 181, 45))
    plt.xlim([-180, 180])
    plt.xlabel("Longitude [degrees]", labelpad=8)
    # yticks
    plt.yticks(range(-90, 91, 30), range(-90, 91, 30))
    plt.ylim([-90, 90])
    plt.ylabel("Latitude [degrees]")

#
#   create_figure_slice
#
def create_figure_slice(figid, aspect, yy, yl, yr, slice, title, yvisible, cmin, cmax, levels, v3d):
    # create figure
    print "# Creating slice ..."
    
    # data
    vv = v3d[0,yl:yr,slice,:]
    # depths, heights, mids
    vzz = np.nanmax(v3d[1, :, slice, :], axis = 0)
    vdz = np.nanmax(v3d[2, :, slice, :], axis = 0)
    vdzz = vzz - 0.5*vdz
#    print vzz, vdz, vdzz

    # plot slice
    plt.figure(figid)
    # colormap
    cmap = plt.cm.bone_r
    # contour fill
    p1 = plt.contourf(yy[yl:yr], vdzz, vv.T, cmap=cmap, levels=levels, origin="upper")
    plt.clim(cmin, cmax)
    # x
    plt.xticks([-60, -30, 0, 30, 60], [-60, -30, 0, 30, 60])
    plt.xlabel("Latitude [degrees]", labelpad=8)
    # y
    plt.ylim([0, 5200])
    plt.yticks(vzz, ("50", "", "", "360", "", "790", "1080", "1420", "1810", "2250", "2740", "3280", "3870", "4510", "5200"))
    plt.yticks(visible=yvisible)
    if yvisible: plt.ylabel("Depth [m]")
    # title
    t1 = plt.title(title)
    t1.set_y(1.02)
    # contour lines
    p2 = plt.contour(yy[yl:yr], vdzz, vv.T, levels=p1.levels, linewidths = (1,), colors="k")#, hold="on")
    # axes transform
    plt.gca().invert_yaxis()
    plt.gca().set_aspect(aspect)
    
    # colorbar
    cbar = plt.colorbar(p1, format = '%.1e', fraction = 0.027, pad = 0.03)

#
#   create_figures
#
def create_figures(m3dprefix, modeldir, tracer, v3d):
	# create figures
	print "# Creating figures ..."

	# longitude
	dx = 2.8125
	xx = np.concatenate([np.arange(-180 + 0.5 * dx, 0, dx), np.arange(0 + 0.5 * dx, 180, dx)])
	# latitude
	dy = 2.8125
	yy = np.arange(-90 + 0.5 * dy, 90, dy)
# 	print xx
# 	print yy
	# min, max of v3d data, just info
	vmin = np.nanmin(v3d[0,:,:,:])
	vmax = np.nanmax(v3d[0,:,:,:])
	vavg = np.nanmean(v3d[0,:,:,:])
	vstd = np.nanstd(v3d[0,:,:,:])
	print("min: %e, max: %e, avg: %e, std: %e" % (vmin, vmax, vavg, vstd))

	# levels
	levels = np.arange(vmin, vmax, vstd)
	print levels

	# color range
	cmin, cmax = [vmin, vmax]
#    print cmin, cmax

    #
    #   slices
    #

	# range of latitude
	yl, yr = [6, 58]
	# aspect of axes
	aspect = (yy[yr] - yy[yl] + 1) / 5200.0 / 1.8
#    print yl, yr
#    print aspect

	slices = [73, 117, 32]
	title = "Pacific"
	yvisible = True
	create_figure_slice(1, aspect, yy, yl, yr, slices[0], title, yvisible, cmin, cmax, levels, v3d)
    # write to file
	plt.savefig("figures/slice-%s-%s-diff-%s.pdf" % (modeldir, tracer, title), bbox_inches = "tight")

	title = "Atlantic"
	yvisible = True
	create_figure_slice(2, aspect, yy, yl, yr, slices[1], title, yvisible, cmin, cmax, levels, v3d)
	# write to file
	plt.savefig("figures/slice-%s-%s-diff-%s.pdf" % (modeldir, tracer, title), bbox_inches = "tight")

	title = "Indian"
	yvisible = True
	create_figure_slice(3, aspect, yy, yl, yr, slices[2], title, yvisible, cmin, cmax, levels, v3d)
	# write to file
	plt.savefig("figures/slice-%s-%s-diff-%s.pdf" % (modeldir, tracer, title), bbox_inches = "tight")

    #
    #   surface
    #
	aspect = 1.0
	create_figure_surface(4, aspect, xx, yy, cmin, cmax, levels, slices, v3d)
    # write to file
	plt.savefig("figures/surface-%s-diff-%s.pdf" % (modeldir, tracer), bbox_inches = "tight")

#
#   main
#
if __name__ == "__main__":
	# no arguments?
	if len(sys.argv) <= 2:
		# print usage and exit with code 1
		print("usage: %s [MODELNAME...] [TRACER...] [PREFIX...]" % sys.argv[0])
		sys.exit(1)
	# model directory does not exist?
	modeldir = sys.argv[1]
	if not (os.path.exists(modeldir) and os.path.isdir(modeldir)):
		# print error message and exit with code 2
		print("### ERROR ### '%s' does not exist or is no directory." % modeldir)
		sys.exit(2)
	# strip trailing slash if any
	modeldir = os.path.normpath(modeldir)
	# debug
	print("using %s model" % modeldir)
	# set prefix for data
	m3dprefix = os.path.expanduser('~/.metos3d')
	# set tracer
	tracer = sys.argv[2]
	# check if we have a *given* prefix
	if len(sys.argv) == 4:
		prefix = sys.argv[3]
	else:
		prefix = "newton."
	# debug
	print("using prefix ... '%s'" % prefix)
	# read data
	v3d = read_data(m3dprefix, modeldir, prefix, tracer)
	# create figure
	create_figures(m3dprefix, modeldir, tracer, v3d)
		
	

