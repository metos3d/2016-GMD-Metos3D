#!/usr/bin/env python

import os
import sys
import re
import numpy as np

#
#	read_monitor_data
#
def read_monitor_data(filepath):
	# debug
# 	print 'read_monitor_data ...', filepath
    # read file, parse and store
	splist = []
	f = open(filepath, 'r')
	for line in f:
		spmatch = re.search(r"^\s+(\d+.\d+)s \d+ Spinup Function norm \d+.\d+e[+-]\d+ \d+.\d+e[+-]\d+", line)
		if spmatch:
# 			print line,
			splist.append(float(spmatch.groups()[0]))
	f.close()
	# return values
	return np.array(splist)

#
#	create table
#
def create_table(models, data):
	# debug
	print("create table ...")
	for i in range(0, len(models)):
		# debug
# 		print models[i], data[i]
		# model
		print("%-15s" % models[i]),
		# min, max, avg, std
		min = data[i][0]
		max = data[i][1]
		avg = data[i][2]
		std = data[i][3]
		print(r''' & '''),
		print("%8.2f" % min),
		print(r'''\unit{s} '''),
		print(r''' & '''),
		print("%8.2f" % max),
		print(r'''\unit{s} '''),
		print(r''' & '''),
		print("%8.2f" % avg),
		print(r'''\unit{s} '''),
		print(r''' & '''),
		print("%8.2f" % std),
		print(r''' & '''),
		print("%8.2f" % (min/(i+1))),
		print(r'''\unit{s} '''),
		# row end
		print(r''' \\''')

# 		print(r'''
# %tab:norms
# \begin{table*}
# \setlength{\tabcolsep}{10pt}
# \caption{\label{tab:norms}
# Difference in Euclidean norm ...}
# \vspace{0.3cm}
# \small
# \centering
# \begin{tabular}{l S S}
# \tophline
# Model			&	{$\| \mathbf{y}_{\mathrm{SP}} - \mathbf{y}_{\mathrm{NK}} \|_2$}	& 	{$\| \mathbf{y}_{\mathrm{SP}} - \mathbf{y}_{\mathrm{NK}} \|_{2, \Omega}$}	\\
# \middlehline''')
# 	# go through models
# 	for i in range(0, len(modellist)):
# 		# combine volumes and tracer vectors
# 		vols = np.tile(vol, len(tracerlist[i]))
# 		tracer = np.concatenate(vlist[i])
# 		# debug
# 		print("%-15s" % modellist[i]),
# 		print(r'''& '''),
# # 		print(r'''& \num{'''),
# # 		print("vector lengths: %d, %d" % (len(vols), len(tracer)))
# 		# compute 2-norm
# 		norm = np.linalg.norm(tracer)
# 		# debug
# 		print("%.3e" % norm),
# 		print(r''' & '''),
# # 		print(r'''} & \num{'''),
# 		# compute weighted norm
# 		normw = np.linalg.norm(tracer * np.sqrt(vols))
# 		# debug
# 		print("%.3e" % normw),
# 		print(r''' \\''')
# # 		print(r'''} \\''')
# 	print(r'''\bottomhline
# \end{tabular}
# \end{table*}
# ''')

#
#   read_all_data
#
def read_all_data(modeldirs):
	data = []
	for dir in modeldirs:
		# construct file path
		filepath = "%s/work/profile.%s.out" % (dir, dir)
		# read monitor data
		monitordata = read_monitor_data(filepath)
		# debug
# 		print monitordata
		# compute min, max, avg, std
		min = np.min(monitordata)
		max = np.max(monitordata)
		avg = np.mean(monitordata)
		std = np.std(monitordata)
		# debug
# 		print "%.2f s" % min, "%.2f s" % max, "%.2f s" % avg, "%.2f" % std
		# append to data
		data.append([min, max, avg, std])
	# return values
	return data

#
#   main
#
if __name__ == "__main__":
	# no arguments?
	if not (len(sys.argv) == 1):
		# print usage and exit with code 1
		print("usage: ./create-table.py")
		sys.exit(1)
	# model directories
	modeldirs = ["N", "N-DOP", "NP-DOP", "NPZ-DOP", "NPZD-DOP"]
	# debug
	print("using %s models" % modeldirs)
	# read all data
	data = read_all_data(modeldirs)
	# debug
# 	print data
	# create table
	create_table(modeldirs, data)








###
###	DUMP
###

	# create profile plot
# 	create_profile_figure(modeldir)

#
#   read_profile_data
#
# def read_profile_data(filepath):
# 	# debug
# 	print 'read profile data ...', filepath
# 	# read file, parse and store
# 	lines = []
# 	starthere = 0
# 	f = open(filepath, "r")
# 	for line in f:
# 		# look for starting point
# 		linematch = re.search(r"--- Event Stage 0: Main Stage", line)
# 		if linematch:
# 			starthere = 1
# 		if starthere == 1:
# 			linematch = re.search(r"------------------------------------------------------------------------------------------------------------------------", line)
# 			if linematch:
# 				break
# 			lines.append(line)
# 	f.close()
# # 	print lines
#
# 	# post process
# 	data = []
# 	total = 0.0
# 	for line in lines:
# 		datamatch = re.search(r"^(\S+)\s+\d+\s\S+\s(\d+.\d+e[+-]\d+)", line)
# 		if datamatch:
# 			groups = datamatch.groups()
# 			if groups[0] in ["BGCStep", "MatCopy", "MatScale", "MatAXPY", "MatMult"]:
# 				data.append(datamatch.groups())
# 			if groups[0] in ["TimeStepPhi"]:
# 				total = float(groups[1])
# 	# return values
# 	return total, data

		# read profile data
# 		profiledata = read_profile_data(filepath)

# 		print profiledata
