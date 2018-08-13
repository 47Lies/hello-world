from lxml import etree
import os
import re
configfile: "config.yaml"

MQ = "/bioinfo/users/ltaing/ProteoGenomics/PROTEOME_GENERATOR/SOFTS/proteomegenerator/MaxQuant/bin/MaxQuantCmd.exe"
MaxQuantWorkingDirectory=config["MaxQuant"]["WorkingDirectory"]

if "PBS_JOBID" in os.environ:
	ScratchMQPAR_Name=MaxQuantWorkingDirectory+"/mqpar"+os.environ["PBS_JOBID"]+".xml"
else:
	ScratchMQPAR_Name=MaxQuantWorkingDirectory+"/mqpar.xml"

def CreateOrClean(DirectoryName):
	if os.path.isdir(DirectoryName):
		print("Directory",DirectoryName,"exist\nProceed to cleaning")
		CMD="rm --recursive --verbose --force "+DirectoryName+"/*"
		os.system(CMD)
	else :
		os.makedirs(DirectoryName)

rule all:
	input: "Toto.txt"

rule MakeScratchyMQPAR_XML:
	input: 
		template="mqparTemplateSILAC.xml"
	output:
		ScratchMQPAR_XML=ScratchMQPAR_Name
	run:
		UTF8_XML_PARSER=etree.XMLParser(encoding='UTF-8',remove_blank_text=True)
		tree=etree.parse(input.template,UTF8_XML_PARSER)
		
		Fasta = tree.xpath("/MaxQuantParams/fastaFiles/FastaFileInfo/fastaFilePath")[0]
		Fasta.text=config["MaxQuant"]["FASTA_DATABASE"]
		
		JobFSF=MaxQuantWorkingDirectory+"/FixedSearchFolder"
		JobFCF=MaxQuantWorkingDirectory+"/fixedCombinedFolder"
		FSF = tree.xpath("/MaxQuantParams/fixedSearchFolder")[0]
		FCF = tree.xpath("/MaxQuantParams/fixedCombinedFolder")[0]
		FSF.text=JobFSF
		FCF.text=JobFCF
		CreateOrClean(JobFSF)
		CreateOrClean(JobFCF)
		
		if "PBS_JOBID" in os.environ:
			TempDir="/local/scratch/"+os.environ["PBS_JOBID"]
		else:
			TempDir=MaxQuantWorkingDirectory+"/TempFolder"
			CreateOrClean(TempDir)
		TF =  tree.xpath("/MaxQuantParams/tempFolder")[0]
		TF.text=TempDir
		
		N =  tree.xpath("/MaxQuantParams/numThreads")[0]
		N.text=config["MaxQuant"]["Threads"]
		
		RawFiles = tree.xpath("/MaxQuantParams/filePaths")[0]
		RawFiles.clear()
		Fraction = tree.xpath("/MaxQuantParams/fractions")[0]
		Fraction.clear()
		Experiments = tree.xpath("/MaxQuantParams/experiments")[0]
		Experiments.clear()
		PTMS = tree.xpath("/MaxQuantParams/ptms")[0]
		PTMS.clear()
		PGI = tree.xpath("/MaxQuantParams/paramGroupIndices")[0]
		PGI.clear()
		NumFile=0
		NumFraction=1
		CurrentExperiment=""
		Pattern=re.compile(".*raw$")
		for root, dirs, files in os.walk(config["RAW_GLOBAL_GIR"]):
			dirs.sort()
			files.sort()
			for file in files:
				if not Pattern.match(file):
					continue
				else:
					p=os.path.join(root,file)
					line=os.path.abspath(p)
					line=line.rstrip()
					RmIndex=line
					RmIndex=RmIndex.replace("raw","index")
					RmIndex="rm --force"+RmIndex
					RmFolder=line
					RmFolder=RmFolder.replace(".raw","")
					RmFolder="rm -R --force"+RmFolder
					DirsFile=line.split('/')
					NameIdx=len(DirsFile)
					Name=DirsFile[NameIdx-1]
					DirName=DirsFile[NameIdx-2]
					if DirName != CurrentExperiment:
						NumFraction=1
						CurrentExperiment=DirName
					else:
						NumFraction=NumFraction+1
					ZefilePath=etree.Element("string")
					ZefilePath.text=line
					RawFiles.append(ZefilePath)
					Zefraction=etree.Element("short")
					Zefraction.text=str(NumFraction)
					Fraction.append(Zefraction)
					ZeExperiment=etree.Element("string")
					ZeExperiment.text=CurrentExperiment
					Experiments.append(ZeExperiment)
					ZePTMS=etree.Element("boolean")
					ZePTMS.text="False"
					PTMS.append(ZePTMS)
					
					ZePGI=etree.Element("int")
					ZePGI.text="0"
					PGI.append(ZePGI)
					
					NumFile=NumFile+1
		tree.write(output.ScratchMQPAR_XML,xml_declaration=True,encoding='UTF-8',pretty_print=True)

rule RunMaxQuant:
	input: 
		ScratchMQPAR_Name
	output: 
		"Toto.txt"
	threads: 24
	shell:
		"mono {MQ} {ScratchMQPAR_Name} > Toto.txt"
