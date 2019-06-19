import FileDesc  
  
testfiles = [
	["test.txt","",""],
	["test2.txt","",""]
	]
	
import subprocess
import os
import getpass

#subprocess.call(["dir"],shell=True)

#get password
password = getpass.getpass()
#print(password)
username = "gustajo@chalmers.se"

			
for elem in FileDesc.files:
	filename = elem[0]
	fileId = elem[1]
	first = True
	while True:
		exists = os.path.isfile(filename)
		if exists:
			sz = os.path.getsize(filename)
			if sz==0:
				os.remove(filename)
				exists = False
		filenameCip = filename + ".cip"
		filenameGpg = filename + ".gpg"
		existsGpg = os.path.isfile(filenameGpg)
		existsCip = os.path.isfile(filenameCip)
		if existsCip:
			sz = os.path.getsize(filenameCip)
			if sz==0:
				os.remove(filenameCip)
				existsCip = False
		if existsGpg:
			sz = os.path.getsize(filenameGpg)
			if sz==0:
				os.remove(filenameGpg)
				existsGpg = False
		tmpfn = filename + ".cip.egastream"
		#always delete the egastream file
		existsTmp = os.path.isfile(tmpfn)
		if existsTmp:
			os.remove(tmpfn)
		tmpfn2 = filename + ".gpg.egastream"
		#always delete the egastream file
		existsTmp2 = os.path.isfile(tmpfn2)
		if existsTmp2:
			os.remove(tmpfn2)
		if exists:
			if first:
				print("Skipping file since it already exists: ", elem[0])
			break
		else:
			print("Trying to download file: ", elem[0])
			
			if not (existsCip or existsGpg):
				#create request
				#    java -jar EgaDemoClient.jar -p gustajo@chalmers.se pwd -rf EGAF00000888764 -re ger -label request_EGAF00000888764
				requestId = "request_" + fileId
				subprocess.call(["java", "-jar", "EgaDemoClient.jar", "-p", username, password,"-rf",fileId,"-re", "ger", "-label",requestId],shell=True)
				
				#download
				#java -jar EgaDemoClient.jar -p gustajo@chalmers.se pwd -dr request_EGAF00000888764 -nt 4
				subprocess.call(["java", "-jar", "EgaDemoClient.jar", "-p", username, password,"-dr",requestId, "-nt", "4"],shell=True)
			else:
				print("Skipping download since .cip/.gpg file already exists: ", elem[0],".cip")

			#Check if cip or gpg
			existsGpg = os.path.isfile(filenameGpg)
			#existsCip = os.path.isfile(filenameCip)
			ext = ".cip";
			if existsGpg:
				ext = ".gpg"
			#decrypt
			#java -jar EgaDemoClient.jar -p email password -dc EGAR00001356280_150616_SN546_0278_B_C794GACXX_AGTTCC_5_read1.fastq.gz.cip EGAR00001356280_150616_SN546_0278_B_C794GACXX_AGTTCC_5_read1.fastq.gz -dck ger
			subprocess.call(["java", "-jar", "EgaDemoClient.jar", "-p", username, password,"-dc",filename + ext, "-dck", "ger"],shell=True)
			
			
		first = False

