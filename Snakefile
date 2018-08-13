rule MakeScratchyMQPAR_XML:
	input: 
		template="mqparTemplateSILAC.xml"
	output:
		ScratchMQPAR_XML="Momo.xml"
	run:
		from lxml import etree
		import os
		import re
		UTF8_XML_PARSER=etree.XMLParser(encoding='UTF-8',remove_blank_text=True)
		tree=etree.parse(input.template,UTF8_XML_PARSER)
		Fasta = tree.xpath("/MaxQuantParams/fastaFiles/FastaFileInfo/fastaFilePath")[0]
		Fasta.text="Paris a le blues"
		FSF =  tree.xpath("/MaxQuantParams/fixedSearchFolder")[0]
		TF =  tree.xpath("/MaxQuantParams/tempFolder")[0]
		N =  tree.xpath("/MaxQuantParams/numThreads")[0]
		FCF =  tree.xpath("/MaxQuantParams/fixedCombinedFolder")[0]
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
		for root, dirs, files in os.walk('/data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/'):
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
	

