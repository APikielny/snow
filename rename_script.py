import os
import argparse
import sys
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory", help="REQUIRED. Path to directory to rename. Relative to the path containing the script")
parser.add_argument("-b", "--basename", help="REQUIRED. Base name of the files to rename. Used so if there are files in the directory that should not be touched, they are not. Example: filename - scene1.jpg.001, the basename is 'scene1.")
args = parser.parse_args()

idx = [0, 2, 1]
base_path = os.getcwd()

for filename in os.listdir(args.directory):
	name = filename.split(os.extsep)
	if len(name) == 3 and name[0] == args.basename:
		newName = np.array(name)
		newName = newName[idx]
		newName = newName.tolist()
		filenameNew = '.'.join(newName)
		relativePathOld = args.directory  + "/" + filename
		relativePathNew = args.directory  + "/" + filenameNew
		os.rename(relativePathOld, relativePathNew)

	

