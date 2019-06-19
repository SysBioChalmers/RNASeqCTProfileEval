import FileDesc	
import os
import subprocess

import shutil

			
for project in FileDesc.kallistoInfo:
	projDir = project[0] # same as project id
	print("project: ", projDir)
	projDirExists = os.path.isdir(projDir)
	if not projDirExists:
		print("Failed: Project dir does not exist: ", projDir)
		break;
	os.chdir(projDir);

	for sample in project[1]:
		file1 = sample[0]
		file2 = sample[1]
		outputDir = sample[2]
		print("sample: ", outputDir)
		outputDirExists = os.path.isdir(outputDir)
		if not outputDirExists: #otherwise we assume this file has been processed
			os.mkdir(outputDir)
			#now call kallisto
			subprocess.call(["..\kallisto.exe", "quant", "-i", "..\\transcripts.idx", "-o", outputDir, "-b", "1", file1, file2],shell=True)
		else:
			print("skipped running sample, the output dir already exists, we assume it has already been run")
		
	os.chdir("..");
	
#the line below needs to be run once to create the transcripts.idx file
#kallisto index -i transcripts.idx Homo_sapiens.GRCh38.cdna.all.fa.gz

#typical case for running kallisto
#kallisto index -i transcripts.idx Homo_sapiens.GRCh38.cdna.all.fa.gz