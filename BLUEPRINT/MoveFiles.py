import FileDesc	
import os

import shutil

			
for elem in FileDesc.files:
	filename = elem[0]
	fileId = elem[1]
	dir = elem[2] # same as project id
	#print(filename," - ", dir)
	dirExists = os.path.isdir(dir)
	if not dirExists:
		os.mkdir(dir)
	fileExists = os.path.isfile(filename);
	if fileExists:
		shutil.move(filename, dir + "/" + filename)
	
