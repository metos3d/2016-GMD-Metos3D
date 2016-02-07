#!/usr/bin/env python -W ignore

import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.rc("font",**{"family":"sans-serif"})
matplotlib.rcParams.update({'font.size': 14})
matplotlib.rc("text", usetex=True)
matplotlib.use("PDF")
import matplotlib.pyplot as plt

#
#   figure_error
#
def figure_error(message):
    print "#"
    print "# ERROR:", message
    print "#"
    print "# Quitting ..."
    sys.exit(0)

#
#   figure_read_nesh_fe_data
#
def figure_read_nesh_fe_data():
    # number of tests
    # number of runs per test
    # number of timed model years per run
    ntest = 10
    nrun = 16*16
    nyear = 3
    # create array
    data = np.empty(shape=(nrun*nyear, ntest), dtype=float)
    # read nesh-fe data
    print "# Reading nesh-fe files ..."
    for itest in range(0, ntest):
        filename = "work/metos3d.nesh-fe.out.MITgcm-PO4-DOP-" + str(itest) + ".txt"
        print "# ", filename
        # read, parse
        f = open(filename, "r")
        # data1
        fiter = 0
        for line in f:
            matches = re.search(r"^\s+(\d+.\d+)s (\d+) Spinup Function norm (\d+.\d+e[+-]\d+)", line)
            if matches:
#                print line,
                [time, iter, norm] = matches.groups()
                data[fiter, itest] = time
                fiter = fiter + 1
        f.close()
        if not fiter == nrun*nyear:
            figure_error("Mismatch: fiter = " + str(fiter) + ", nrun * nyear = " + str(nrun*nyear) )
    # post process
    mindata = np.empty(shape=(nrun), dtype=float)
    for irun in range(0, nrun):
        mindata[irun] = np.amin(data[irun*nyear:(irun+1)*nyear, :])
    # seq
    data = np.empty(shape=(nyear), dtype=float)
    filename = "work/metos3d.nesh-fe.out.MITgcm-PO4-DOP-seq.txt"
    f = open(filename, "r")
    dataiter = 0
    print "# ", filename
    for line in f:
        matches = re.search(r"^\s+(\d+.\d+)s (\d+) Spinup Function norm (\d+.\d+e[+-]\d+)", line)
        if matches:
            [time, iter, norm] = matches.groups()
            data[dataiter] = time
            dataiter = dataiter + 1
    f.close()
    # post process
    minseqdata = np.amin(data)
    # set t1
    mindata[0] = np.minimum(mindata[0], minseqdata)
	# return data
    return mindata

#
#	figure_read_spk_data
#
def figure_read_spk_data():
	# debug
	print("Reading TMM data ...")
	# read, parse
	f = open("work/spkwall.out", "r")
	# go through lines
	nrun = 5
	irun = 1
	mintime = 1.e+300
	spkdata = []
	for line in f:
		matches = re.search(r"^Wall clock time:\s+(\d+.\d+)", line)
		if matches:
			time = float(matches.groups()[0])
# 			print irun, irun % nrun, mintime, time
			mintime = np.minimum(mintime, time)
			if irun % nrun == 0:
# 				print "append"
				spkdata.append(mintime)
				irun = 1
				mintime = 1.e+300
				continue
			irun = irun + 1
	# close
	f.close()
	# return data
	return spkdata

#
#   figure_create_plot_speedup
#
def figure_create_plot_speedup(nesh_fe_data, spk_data, load_data):
    # create figure
    print "# Creating speedup figure ..."
    # create plot
    plt.figure(1)
    # set range
    nmax = 260
    xrange = np.arange(1, nmax+1)
    # ticks
    plt.xticks(np.arange(0, nmax+1, 20), np.arange(0, nmax+1, 20))
    plt.yticks(np.arange(0, nmax+1, 20), np.arange(0, nmax+1, 20))
    # axis
    plt.axis([0, nmax, 0, nmax])
    # grid
    plt.grid(True, which="major", axis="both", color="0.25", linestyle=":")
    # axis title
    xlab = plt.xlabel("Number of processes", labelpad=10)
    ylab = plt.ylabel("Speedup factor", labelpad=10)
    # ideal speedup, load balance
    p1, = plt.plot(xrange, xrange, "--", lw=1, color="0.0")
    p2, = plt.plot(xrange, xrange * load_data, "-", lw=1, color="0.0")
    # compute and plot speedup
    nesh_fe_data = nesh_fe_data[0] / nesh_fe_data
    spk_data = spk_data[0] / spk_data
    p3, = plt.plot(xrange[0:256], nesh_fe_data, "-", lw=1, color="0.5")
    p4, = plt.plot(xrange[0:192], spk_data, "-", lw=1, color="0.7")
    
    # legend
    leg1 = r"Perfect"
    leg2 = r"Theoretical limit of Metos3D load balancing"
    leg3 = r"Metos3D on Intel$^{\mbox{\textregistered}}$ Sandy Bridge EP hardware"
    leg4 = r"TMM on Intel$^{\mbox{\textregistered}}$ Sandy Bridge EP hardware"
    plt.legend([p1,p2,p3,p4],[leg1,leg2,leg3,leg4], loc=2, numpoints=1, fontsize = 12.0)
    # save
    plt.figure(1)
    plt.savefig("figures/speedup.pdf", bbox_inches="tight")

#
#   figure_create_plot_efficiency
#
def figure_create_plot_efficiency(nesh_fe_data, spk_data, load_data):
	# create figure
	print "# Creating efficiency figure ..."
	# create plot
	plt.figure(2)
	# range
	nmax = 260
	neffmax = 140
	xrange = np.arange(1, nmax+1)
	# axis
	plt.axis([0, nmax, 0, neffmax])
	# ticks
	plt.xticks(np.arange(0, nmax+1, 20), np.arange(0, nmax+1, 20))
	plt.yticks(np.arange(0, neffmax+1, 10), np.arange(0, neffmax+1, 10))
	# grid
	plt.grid(True, which="major", axis="both", color="0.25", linestyle=":")
	# axis title
	xlab = plt.xlabel("Number of processes", labelpad=10)
	ylab = plt.ylabel(r"Efficiency [\%]", labelpad=10)
	# ideal and theoretical efficiency
	p1, = plt.plot(xrange, np.repeat(100.0, nmax), "--", lw=1, color="0.0")
	p2, = plt.plot(xrange, load_data * 100.0, "-", lw=1, color="0.0")
	# compute and plot efficiency
	nesh_fe_data = 100.0 * nesh_fe_data[0] / nesh_fe_data / xrange[0:256]
	spk_data = 100.0 * spk_data[0] / spk_data / xrange[0:192]
	# plot
	p3, = plt.plot(xrange[0:256], nesh_fe_data, "-", lw=1, color="0.5")
	p4, = plt.plot(xrange[0:192], spk_data, "-", lw=1, color="0.7")
	# save
	plt.figure(2)
	plt.savefig("figures/efficiency.pdf", bbox_inches="tight")

#
#   figure_compute_load
#
def figure_compute_load(m3dprefix):
    print "# Computing load balancing ..."
    import numpy
    # read profile data
    # read starts
    filename = "%s/data/data/TMM/2.8/Geometry/gStartIndices.bin" % m3dprefix
    print "#", filename
    f = open(filename, "rb")
    nprof = numpy.fromfile(f, dtype=">i4", count=1)
    starts = numpy.fromfile(f, dtype=">i4", count=nprof)
    f.close()
    # read ends
    filename = "%s/data/data/TMM/2.8/Geometry/gEndIndices.bin" % m3dprefix
    print "#", filename
    f = open(filename, "rb")
    nprof = numpy.fromfile(f, dtype=">i4", count=1)
    ends = numpy.fromfile(f, dtype=">i4", count=nprof)
    f.close()
    # compute lengths
    lengths = ends-starts+1
    veclength = numpy.sum(lengths)
    # compute theoretical load balancing
    nprocs = 260
    theory = numpy.zeros(shape=(nprocs), dtype=float)
    # prepare index computation
    idxwork = (starts - 1 + 0.5 * lengths)/veclength
    # compute load for all procs from 1 to nprocs
    for iproc in range(1, nprocs+1):
        # optimal vectorlength
        optlength = veclength/float(iproc)
        # storage for all lengths
        profbucket = numpy.zeros(shape=(iproc), dtype=">i4")
        # indices for iproc
        idx = numpy.floor(idxwork*iproc)
        # compute balanced load for iproc
        for iprof in range(0, nprof):
            intidx = int(idx[iprof])
            profbucket[intidx] = profbucket[intidx] + lengths[iprof]
        # take worst case
        theory[iproc-1] = veclength / float(numpy.amax(profbucket)*iproc)
    return theory

#
#   figure_main
#
if __name__ == "__main__":
    # read data
	nesh_fe_data = figure_read_nesh_fe_data()
# 	print nesh_fe_data
	spk_data = figure_read_spk_data()
# 	print spk_data
	# set prefix for data
	m3dprefix = os.path.expanduser('~/.metos3d')
	# compute load balancing
	load_data = figure_compute_load(m3dprefix)
    # create figures
	figure_create_plot_speedup(nesh_fe_data, spk_data, load_data)
	figure_create_plot_efficiency(nesh_fe_data, spk_data, load_data)


