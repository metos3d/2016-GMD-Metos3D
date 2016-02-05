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

#
#	read_nk_convergence_data
#
def read_nk_convergence_data(nmax, filepath):
	# debug
	print 'read Newton-Krylov convergence data ...'
	# read file, parse and store
	kspiter = 0
	kspidx = []
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
			kspidx.append(kspiter)
			ksplist.append(kspmatch.groups())
			kspiter = kspiter + 1
			if kspiter >= nmax:
				break
	f.close()
	# debug	
	print("read in %d data points ..." % kspiter)	
	# retrieve values
	kspval = map(lambda x: x[1], ksplist)
	snesval = map(lambda x: x[1], sneslist)
	# return values
	return [kspidx, kspval, snesidx, snesval]

#
#	read_sp_convergence_data
#
def read_sp_convergence_data(nmax, filepath):
	# debug
	print 'read spin-up convergence data ...'
    # read file, parse and store
	spiter = 0
	spidx = []
	splist = []
	f = open(filepath, 'r')
	for line in f:
		spmatch = re.search(r"^\s+\d+.\d+s (\d+) Spinup Function norm (\d+.\d+e[+-]\d+) (\d+.\d+e[+-]\d+)", line)
		if spmatch:
# 			print line,
			spidx.append(spiter)
			splist.append(spmatch.groups())
			spiter = spiter + 1
			if spiter >= nmax:
				break
			else:
				continue
			break
	f.close()
	# debug	
	print("read in %d data points ..." % spiter)
    # retrieve values
	spval = map(lambda x: x[1], splist)
	spvalw = map(lambda x: x[2], splist)
	# return values
	return [spidx, spval, spvalw]

#
#	create_convergence_figure
#
def create_convergence_figure(modeldir, prefix):
	# debug
	print "create convergence figure ..."
	# set max data points
	nmax = 10000
	# read convergence data
	(kspidx, kspval, snesidx, snesval) = read_nk_convergence_data(nmax, "%s/work/%s%s.out" % (modeldir, prefix, modeldir))
	# spinup
	(spidx, spval, spvalw) = read_sp_convergence_data(nmax, "%s/work/spinup.%s.out" % (modeldir, modeldir))	
	
	# set dimensions
	# plot
	xrange = range(0, nmax)
	# subplot
	nsplit = 1000
	ax1l, ax1r = (0, nsplit)
	ax2l, ax2r = (nsplit, nmax)

	# subplots
	figid = 1
	f = plt.figure(figid)
	gs = gsp.GridSpec(1, 2, width_ratios=[2,1])
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1])
	
	# spinup
	p1, = ax1.semilogy(spidx[ax1l:ax1r], spval[ax1l:ax1r], "k-", linewidth = 1.0)
	ax2.semilogy(spidx[ax2l:ax2r], spval[ax2l:ax2r], "k-", linewidth = 1.0)

	# KSP
	p2, = ax1.semilogy(kspidx[ax1l:ax1r], kspval[ax1l:ax1r], "k-", linewidth = 1.0)
	ax2.semilogy(kspidx[ax2l:ax2r], kspval[ax2l:ax2r], "k-", linewidth = 1.0)

	# SNES
	if max(snesidx) >= nsplit:
		# get index lists first NK
		idxlow = np.where(np.array(snesidx)<nsplit)[0]
		idxhigh = np.where(np.array(snesidx)>=nsplit)[0]
		# set arrays first NK
		snesidxlow = [snesidx[i] for i in idxlow]
		snesidxhigh = [snesidx[i] for i in idxhigh]
		snesvallow = [snesval[i] for i in idxlow]
		snesvalhigh = [snesval[i] for i in idxhigh]
		# plot
		ax1.semilogy(snesidxlow, snesvallow, marker = "D", ms = 4.0, mfc = "None", mec = "k", mew = 1.0, linewidth = 0)
		ax2.semilogy(snesidxhigh, snesvalhigh, marker = "D", ms = 4.0, mfc = "None", mec = "k", mew = 1.0, linewidth = 0)
	else:
		ax1.semilogy(snesidx, snesval, marker = "D", ms = 4.0, mfc = "None", mec = "k", mew = 1.0, linewidth = 0)

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
	ax1.set_ylim([10e-7, 10e+2])
	ytext = [r"10\textsuperscript{%d}" % i for i in range(-7, 4)]
	ax1.set_ylabel("Norm [$\mathrm{m\,mol\,P / m^3}$]")
	ax1.set_yticklabels(ytext)

	ax2.set_yticklabels([])
	ax2.set_ylim([10e-7, 10e+2])

	# grid
	ax1.grid(True, which='major', axis='both', color='0.25', linestyle=':')
	ax2.grid(True, which='major', axis='both', color='0.25', linestyle=':')
	ax2.grid(True, which='minor', axis='x', color='0.25', linestyle=':')
	# legend
	leg1 = r"Spin-up"
	leg2 = r"Newton-Krylov"
	l1 = plt.figlegend([p1, p2], [leg1, leg2], 1, numpoints = 3, bbox_to_anchor = (0.862, 0.899), prop={'size':15})
	ll1 = l1.get_lines()[1]
	ll1.set_marker("D")
	ll1.set_ms(4.0)
	ll1.set_mfc("None")
	ll1.set_mew(1.0)
	# adjust subplots
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.00, hspace=None)
	
	# write to file
	plt.savefig("figures/convergence.%s.pdf" % modeldir, bbox_inches="tight")

#
#   main
#
if __name__ == "__main__":
	# no arguments?
	if len(sys.argv) == 1:
		# print usage and exit with code 1
		print("usage: %s [MODELNAME...] [PREFIX...]" % sys.argv[0])
		sys.exit(1)
	# model directory does not exist?
	modeldir = sys.argv[1]
	if not (os.path.exists(modeldir) and os.path.isdir(modeldir)):
		# print error message and exit with code 2
		print("### ERROR ### '%s' does not exist or is no directory." % modeldir)
		sys.exit(2)
	# strip trailing slash if any
	modeldir = os.path.normpath(modeldir)
	# check if we have a *given* prefix
	if len(sys.argv) == 3:
		prefix = sys.argv[2]
	else:
		prefix = "newton."
	# debug
	print("using '%s' model" % modeldir)
	print("using prefix ... '%s'" % prefix)
	# create convergence plot
	create_convergence_figure(modeldir, prefix)

