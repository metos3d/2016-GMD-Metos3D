#!/usr/bin/env python

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

#
#   read_profile_data
#
def read_profile_data(filepath):
	# debug
	print 'read profile data ...'
	# read file, parse and store
	lines = []
	starthere = 0
	f = open(filepath, "r")
	for line in f:
		# look for starting point
		linematch = re.search(r"--- Event Stage 0: Main Stage", line)
		if linematch:
			starthere = 1
		if starthere == 1:
			linematch = re.search(r"------------------------------------------------------------------------------------------------------------------------", line)
			if linematch:
				break
		lines.append(line)
	f.close()
	# post process
	data = []
	total = 0.0
	for line in lines:
		datamatch = re.search(r"^(\S+)\s+\d+\s\S+\s(\d+.\d+e[+-]\d+)", line)
		if datamatch:
			groups = datamatch.groups()
			if groups[0] in ["BGCStep", "MatCopy", "MatScale", "MatAXPY", "MatMult"]:
				data.append(datamatch.groups())
			if groups[0] in ["TimeStepPhi"]:
				total = float(groups[1])
	# return values
	return total, data

#
#	create_profile_figure
#
def create_profile_figure(modeldir):
	# debug
	print "create profile figure ..."
	# read in profile data
	(total, data) = read_profile_data("%s/work/profile.%s.out" % (modeldir, modeldir))

# 	print total, data
	# sort
	order = ["BGCStep", "MatCopy", "MatScale", "MatAXPY", "MatMult"]
	def new_order(t1, t2):
		return order.index(t1[0]) - order.index(t2[0])
	data = sorted(data, cmp = new_order)
# 	print data
	# sum
	sum = 0.0
	for tuple in data:
		sum = sum + float(tuple[1])
# 	print sum
	# split lists, add other
	dataname = []
	datatime = []
	for tuple in data:
		dataname.append(tuple[0])
		datatime.append(float(tuple[1]))
	dataname.append('Other')
	datatime.append(total-sum)
# 	print dataname, datatime

	# plot   
	fig = plt.figure(1)
	ax = fig.add_subplot(111, aspect = "equal")
	colors = plt.cm.bone(np.linspace(0.2, 1, 6))
	patches, texts, autotexts = ax.pie(datatime, labels = dataname, autopct = "%1.1f\,\%%", colors = colors)
	for text in autotexts:
		text.set_bbox(dict(edgecolor = "k", facecolor = "w"))
	
	# save
	plt.savefig("figures/profile.%s.pdf" % modeldir, bbox_inches = "tight")

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
	# strip trailing slash if any
	modeldir = os.path.normpath(modeldir)
	# debug
	print("using %s model" % modeldir)
	# create profile plot
	create_profile_figure(modeldir)
