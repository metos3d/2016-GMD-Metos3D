#!/usr/bin/env python

import os
import sys
import re
import numpy as np

import matplotlib
matplotlib.rc("font", **{"family" : "sans-serif"})
matplotlib.rcParams.update({'font.size': 14})
# matplotlib.rc("text", usetex = True)
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp

########################################################################
### output
########################################################################

# print_usage
def print_usage():
    print("usage: ./create-figure.py [MODELNAME...]")

########################################################################
### PETSc
########################################################################

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
#     import numpy
    
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
    lsmfull = np.zeros(shape=(nx, ny), dtype=int)
    offset = 0
    for i in range(nx):
        if not nrow[i] == 0.0:
            for j in range(nrow[i]):
                lsmfull[i, colidx[offset]] = int(val[offset])
                offset = offset + 1

    return lsmfull

#
#   read_data
#
def read_data():
#     import numpy as np
    # arrays
    
    # v1d, z, dz, lsm (land sea mask)
    v1d = read_PETSc_vec("exp-01/work/10/PO4.petsc")
    z   = read_PETSc_vec("exp-01/data/TMM/2.8/Forcing/DomainCondition/z.petsc")
    dz  = read_PETSc_vec("exp-01/data/TMM/2.8/Forcing/DomainCondition/dz.petsc")
    lsm = read_PETSc_mat("exp-01/data/TMM/2.8/Geometry/landSeaMask.petsc")
    
    # dims
    nx, ny = lsm.shape
    nz = 15
    # v3d
    v3d = np.zeros(shape=(3, nx, ny, nz), dtype=float)
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
    plt.clabel(p2, fmt = "%2.1f", colors = "k", fontsize = 14)
    # axes transform
    plt.gca().invert_yaxis()
    plt.gca().set_aspect(aspect)
    # write to file
    plt.savefig("solution-slice-" + title, bbox_inches="tight")

#
#   create_figure_surface
#
def create_figure_surface(figid, aspect, xx, yy, cmin, cmax, levels, slices, v3d):
    # basemap
    from mpl_toolkits.basemap import Basemap
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
    p1 = plt.contourf(xx, yy, vv, cmap=cmap, levels=levels, origin="lower")#, hold="on")
    plt.clim(cmin, cmax)
    # contour lines
    p2 = plt.contour(xx, yy, vv, levels=levels, linewidths = (1,), colors="k")#, hold="on")
    plt.clabel(p2, fmt = "%2.1f", colors = "k", fontsize = 14)
    # slices
    s1 = xx[np.mod(slices[0]+64, 128)]
    s2 = xx[np.mod(slices[1]+64, 128)]
    s3 = xx[np.mod(slices[2]+64, 128)]
#    print s1, s2, s3
    plt.vlines([s1, s2, s3], -90, 90, color='k', linestyles='--')
    # set aspect ratio of axes
    plt.gca().set_aspect(aspect)

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

    # write to file
    plt.savefig("solution-surface", bbox_inches="tight")

#
#   create_figure_convergence
#
def create_figure_convergence(figid):
    # basemap
    from mpl_toolkits.basemap import Basemap
    print "# Creating convergence figure ..."

    # max readin
    nmax = 10000
    # read file, parse and store
    kspiter = 0
    ksplist = []
    snesidx = []
    sneslist = []
    f = open("exp-02/work/exp-02.out.txt", "r")
    for line in f:
        # SNES
        snesmatch = re.search(r"^\s+(\d+) SNES Function norm (\d+.\d+e[+-]\d+)", line)
        if snesmatch:
#            print snesmatch.groups()
            snesidx.append(kspiter)
            sneslist.append(snesmatch.groups())
        # KSP
        kspmatch = re.search(r"^\s+(\d+) KSP Residual norm (\d+.\d+e[+-]\d+)", line)
        if kspmatch:
#            print kspmatch.groups()
            ksplist.append(kspmatch.groups())
            kspiter = kspiter + 1
            if kspiter >= nmax:
                break
    f.close()

    # read file, parse and store
    nchunk = 10
    spiter = 0
    splist = []
    for ichunk in range(1, nchunk + 1):
        frmt = str(ichunk).zfill(2)
        filename = "exp-01/work/" + frmt + "/exp-01." + frmt + ".out.txt"
#        print filename
        f = open(filename, "r")
        for line in f:
            spmatch = re.search(r"^\s+\d+.\d+s (\d+) Spinup Function norm (\d+.\d+e[+-]\d+)", line)
            if spmatch:
#                print line,
                splist.append(spmatch.groups())
                spiter = spiter + 1
                if spiter >= nmax:
                    break
        else:
            continue
        break
        f.close()

    # read file, parse and store
    trajiter = 0
    trajlist = []
    f = open("exp-01/work/exp-01.out.txt", "rb")
    for line in f:
#        print line,
        trajmatch = re.search(r"^\d+ \d+ \d+.\d+e[+-]\d+ (\d+.\d+e[+-]\d+)", line)
        if trajmatch:
            trajlist.append(trajmatch.groups())
            trajiter = trajiter + 1
            if trajiter >= nmax:
                break
    f.close()

    # retrieve values
    kspval = map(lambda x: x[1], ksplist)
    snesval = map(lambda x: x[1], sneslist)
    spval = map(lambda x: x[1], splist)
    trajval = map(lambda x: x[0], trajlist)

    # subplots
    f = plt.figure(figid)
    gs = gsp.GridSpec(1, 2, width_ratios=[2,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    # plot
    xrange = range(0, nmax)
    ax1l, ax1r = (0, 1000)
    ax2l, ax2r = (1000, 10000)
    # spinup
    p1, = ax1.semilogy(xrange[ax1l:ax1r], spval[ax1l:ax1r], "k-", linewidth = 1.0)
    ax2.semilogy(xrange[ax2l:ax2r], spval[ax2l:ax2r], "k-", linewidth = 1.0)
    # traj
    p2, = ax1.semilogy(range(100, ax1r+1, 100), trajval[0:10], "k:", marker="o", ms = 4.0, linewidth = 1.0)
    ax2.semilogy(range(1000, ax2r+1, 100), trajval[9:100], "k:", marker="o", ms = 2.0, linewidth = 1.0)
#    p2, = ax1.semilogy(range(100, ax1r+1, 100), trajval[0:10], "k", marker="o", ms = 4.0, linewidth = 1.0)
#    ax2.semilogy(range(1000, ax2r+1, 100), trajval[9:100], "k", marker="o", ms = 2.0, linewidth = 1.0)
    # KSP
    p3, = ax1.semilogy(xrange[ax1l:ax1r], kspval[ax1l:ax1r], "k-", linewidth = 1.0)
    ax2.semilogy(xrange[ax2l:ax2r], kspval[ax2l:ax2r], "k-", linewidth = 1.0)
    # SNES
    ax1.semilogy(snesidx, snesval, color = "k", marker = "D", ms = 4.0, mfc = "None", mec = "k", mew = 1.0, linewidth = None)
    # x
    ax1xmaj = range(0, 1000, 100)
    ax1.set_xticks(ax1xmaj)
    ax1.set_xticklabels(ax1xmaj)
    ax2xmaj = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
    ax2.set_xticks(ax2xmaj)
    ax2.set_xticklabels(ax2xmaj, **{"rotation":55.0, "rotation_mode":"anchor", "x":0.005, "y":0.005, "ha":"right"})

    # common x label
#     f.suptitle("Model years", y = 0.05)
    f.suptitle("Model years", y = 0.035)
#     f.suptitle("Model years")
                    
    # y
    ytext = [r"10\textsuperscript{-5}", r"10\textsuperscript{-4}", r"10\textsuperscript{-3}",
             r"10\textsuperscript{-2}", r"10\textsuperscript{-1}", r"10\textsuperscript{0}",
             r"10\textsuperscript{1}", r"10\textsuperscript{2}"]
    ax1.set_yticklabels(ytext)
    ax1.set_ylabel("Norm [$\mathrm{m\,mol\,P / m^3}$]")
    ax2.set_ylim([10e-6, 10e+1])
    ax2.set_yticklabels([])
    # grid
    ax1.grid(True, which='major', axis='both', color='0.25', linestyle=':')
    ax2.grid(True, which='major', axis='both', color='0.25', linestyle=':')
    ax2.grid(True, which='minor', axis='x', color='0.25', linestyle=':')
    # legend
    leg1 = r"Spin-up, $|| \mathbf{y}_{l} - \mathbf{y}_{l-1} ||_2$"
    leg2 = r"Spin-up, $|| \mathbf{y}_{l} - \mathbf{y}_{l-100} ||_{2,I \times \Omega}$"
    leg3 = r"Newton-Krylov"
#    leg3 = r"Newton-Krylov, $|| \mathbf{y}_{l,j=0} - \mathbf{y}_{l+1,j=0} ||_2$"
    l1 = plt.figlegend([p1, p2, p3], [leg1, leg2, leg3], 1, numpoints = 3, bbox_to_anchor = (0.862, 0.899), prop={'size':15})
    ll1 = l1.get_lines()[2]
    ll1.set_marker("D")
    ll1.set_ms(4.0)
    ll1.set_mfc("None")
    ll1.set_mew(1.0)
    # adjust subplots
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.00, hspace=None)
    
    # write to file
    plt.savefig("solution-convergence", bbox_inches="tight")

#
#   create_figures
#
def create_figures(v3d):
#     import numpy as np
    # create figures
    print "# Creating figures ..."

    # longitude
    dx = 2.8125
    xx = np.concatenate([np.arange(-180 + 0.5 * dx, 0, dx), np.arange(0 + 0.5 * dx, 180, dx)])
    # latitude
    dy = 2.8125
    yy = np.arange(-90 + 0.5 * dy, 90, dy)
#    print xx
#    print yy
#    # min, max of v3d data, just info
#    vmin = np.nanmin(v3d[0,:,:,:])
#    vmax = np.nanmax(v3d[0,:,:,:])
#    print vmin, vmax
    # color range
    cmin, cmax = [0.0, 4.5]
#    print cmin, cmax

    #
    #   slices
    #

    # range of latitude
    yl, yr = [6, 58]
    # aspect of axes
    aspect = (yy[yr] - yy[yl] + 1) / 5200.0 / 1.8
    # levels
    levels = np.arange(0.0, 3.3, 0.4)
#    print yl, yr
#    print aspect
#    print levels

    slices = [73, 117, 32]
#     title = "Pacific"
#     yvisible = True
#     create_figure_slice(1, aspect, yy, yl, yr, slices[0], title, yvisible, cmin, cmax, levels, v3d)
#     title = "Atlantic"
#     yvisible = True
#     create_figure_slice(2, aspect, yy, yl, yr, slices[1], title, yvisible, cmin, cmax, levels, v3d)
#     title = "Indian"
#     yvisible = True
#     create_figure_slice(3, aspect, yy, yl, yr, slices[2], title, yvisible, cmin, cmax, levels, v3d)

    #
    #   surface
    #
    aspect = 1.0
#     create_figure_surface(4, aspect, xx, yy, cmin, cmax, levels, slices, v3d)

    #
    #   convergence
    #
    create_figure_convergence(5)

#
#	read_nk_convergence_data
#
def read_nk_convergence_data(nmax, filepath):
	# debug
	print 'read Newton-Krylov convergence data ...'
	# read file, parse and store
	kspiter = 0
	ksplist = []
	snesidx = []
	sneslist = []
	f = open(filepath, 'r')
	for line in f:
		# SNES
		snesmatch = re.search(r"^\s+(\d+) SNES Function norm (\d+.\d+e[+-]\d+)", line)
		if snesmatch:
# 			print snesmatch.groups()
			snesidx.append(kspiter)
			sneslist.append(snesmatch.groups())
		# KSP
		kspmatch = re.search(r"^\s+(\d+) KSP Residual norm (\d+.\d+e[+-]\d+)", line)
		if kspmatch:
# 			print kspmatch.groups()
			ksplist.append(kspmatch.groups())
			kspiter = kspiter + 1
			if kspiter >= nmax:
				break
	f.close()
	# retrieve values
	kspval = map(lambda x: x[1], ksplist)
	snesval = map(lambda x: x[1], sneslist)
	# return values
	return [kspval, snesidx, snesval]

#
#	read_sp_convergence_data
#
def read_sp_convergence_data(nmax, filepath):
	# debug
	print 'read spin-up convergence data ...'
    # read file, parse and store
	spiter = 0
	splist = []
	f = open(filepath, 'r')
	for line in f:
		spmatch = re.search(r"^\s+\d+.\d+s (\d+) Spinup Function norm (\d+.\d+e[+-]\d+) (\d+.\d+e[+-]\d+)", line)
		if spmatch:
# 			print line,
			splist.append(spmatch.groups())
			spiter = spiter + 1
			if spiter >= nmax:
				break
			else:
				continue
			break
	f.close()
    # retrieve values
	spval = map(lambda x: x[1], splist)
	spvalweight = map(lambda x: x[2], splist)
	# return values
	return [spval, spvalweight]

#
#	create_convergence_plot
#
def create_convergence_plot(figid, nmax, nk_conv_data, sp_conv_data, filepath):
	# debug
	print "create convergence figure ..."
	# retrieve values
	kspval = nk_conv_data[0]
	snesidx = nk_conv_data[1]
	snesval = nk_conv_data[2]
	spval = sp_conv_data[0]
	spvalweight = sp_conv_data[1]
	# subplots
	f = plt.figure(figid)
	gs = gsp.GridSpec(1, 2, width_ratios=[2,1])
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1])
	# plot
	xrange = range(0, nmax)
	ax1l, ax1r = (0, 1000)
	ax2l, ax2r = (1000, nmax)
	# spinup
	p1, = ax1.semilogy(xrange[ax1l:ax1r], spval[ax1l:ax1r], "k-", linewidth = 1.0)
	ax2.semilogy(xrange[ax2l:ax2r], spval[ax2l:ax2r], "k-", linewidth = 1.0)
	p2, = ax1.semilogy(xrange[ax1l:ax1r], spvalweight[ax1l:ax1r], "k-", linewidth = 1.0)
	ax2.semilogy(xrange[ax2l:ax2r], spvalweight[ax2l:ax2r], "k-", linewidth = 1.0)
	# KSP
	p3, = ax1.semilogy(xrange[ax1l:ax1r], kspval[ax1l:ax1r], "k-", linewidth = 1.0)
	ax2.semilogy(xrange[ax2l:ax2r], kspval[ax2l:ax2r], "k-", linewidth = 1.0)
	# SNES
	ax1.semilogy(snesidx, snesval, color = "k", marker = "D", ms = 4.0, mfc = "None", mec = "k", mew = 1.0, linewidth = None)
	# x
	ax1xmaj = range(0, 1000, 100)
	ax1.set_xticks(ax1xmaj)
	ax1.set_xticklabels(ax1xmaj)
	ax2xmaj = range(1000, nmax+1, 1000)
	ax2.set_xticks(ax2xmaj)
	ax2.set_xticklabels(ax2xmaj, **{"rotation":55.0, "rotation_mode":"anchor", "x":0.005, "y":0.005, "ha":"right"})
	# common x label
	f.suptitle("Model years", y = 0.035)
	# y
	ax1.set_ylim([10e-6, 10e+8])
	ax2.set_ylim([10e-6, 10e+8])
	ax2.set_yticklabels([])
	
	
# 	ytext = [r"10\textsuperscript{-5}", r"10\textsuperscript{-4}", r"10\textsuperscript{-3}",
# 			 r"10\textsuperscript{-2}", r"10\textsuperscript{-1}", r"10\textsuperscript{0}",
# 			 r"10\textsuperscript{1}", r"10\textsuperscript{2}"]
# 	ax1.set_yticklabels(ytext)
# 	ax1.set_ylabel("Norm [$\mathrm{m\,mol\,P / m^3}$]")
# 	ax2.set_ylim([10e-6, 10e+1])
# 	ax2.set_yticklabels([])
	# grid
	ax1.grid(True, which='major', axis='both', color='0.25', linestyle=':')
	ax2.grid(True, which='major', axis='both', color='0.25', linestyle=':')
	ax2.grid(True, which='minor', axis='x', color='0.25', linestyle=':')
	# legend
	leg1 = r"Spin-up, $|| \mathbf{y}_{l} - \mathbf{y}_{l-1} ||_2$"
	leg2 = r"Spin-up, $|| \mathbf{y}_{l} - \mathbf{y}_{l-100} ||_{2,I \times \Omega}$"
	leg3 = r"Newton-Krylov"
# 	l1 = plt.figlegend([p1, p2, p3], [leg1, leg2, leg3], 1, numpoints = 3, bbox_to_anchor = (0.862, 0.899), prop={'size':15})
# 	ll1 = l1.get_lines()[2]
# 	ll1.set_marker("D")
# 	ll1.set_ms(4.0)
# 	ll1.set_mfc("None")
# 	ll1.set_mew(1.0)
	# adjust subplots
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.00, hspace=None)
	# write to file
	plt.savefig(filepath, bbox_inches="tight")
	
#
#   main
#
if __name__ == "__main__":
	# no arguments?
	if len(sys.argv) == 1:
		# print usage and exit with code 1
		print("usage: ./create-figure.py [MODELNAME...]")
		sys.exit(1)
	# model directory does not exist?
	modeldir = sys.argv[1]
	if not (os.path.exists(modeldir) and os.path.isdir(modeldir)):
		# print error message and exit with code 2
		print("### ERROR ### '%s' does not exist or is no directory." % modeldir)
		sys.exit(2)


	# set prefix
	m3dprefix = os.path.expanduser('~/.metos3d')
	# set max data points
	nmax = 10000
	# read convergence data
	nk_conv_data = read_nk_convergence_data(nmax, '../work/newton.MITgcm-PO4-DOP.out')
	sp_conv_data = read_sp_convergence_data(nmax, '../work/spinup.MITgcm-PO4-DOP.out')
	# create conv plot
	create_convergence_plot(1, nmax, nk_conv_data, sp_conv_data, 'convergence.MITgcm-PO4-DOP.pdf')
	
# 	print sp_conv_data
	
	
#    z, dz, v3d = read_data()
#     v3d = read_data()
#     create_figures(v3d)

