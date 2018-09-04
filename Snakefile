from lxml import etree
from bs4 import BeautifulSoup
from urllib import request
from ftplib import FTP
import os
import re
os.environ['LANG'] = "en_US.UTF-8"
configfile: "config.yaml"
MQ = "/data/users/ltaing/ProteoGenomics/MaxQuant/Version1.6.2.6/bin/MaxQuantCmd.exe"
MaxQuantWorkingDirectory=config["MaxQuant"]["WorkingDirectory"]
Threads=config["MaxQuant"]["Threads"]

Strawberry=config["Strawberry"]["StrawPath"]

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

DD=BeautifulSoup(request.urlopen('https://www.gencodegenes.org/releases/current.html'),features="lxml")
version_HTML=str(DD.find("h1"))
#'<h1>Release 28 (GRCh38.p12)</h1>' something like that
pre_GencodeVersion=re.sub("^[^ ]+ ","",version_HTML) #discard anything before the first space included starting from the begining
GencodeVersion=re.sub(" .*$","",pre_GencodeVersion) #discard anything after the first space included starting from the end

pre_GenomeVersion=re.sub("^[^(]+\(","",version_HTML) #discard anything before the first parenthesis starting from the begining
GenomeVersion=re.sub("\).*$","",pre_GenomeVersion) # discard anything after the the first parenthesis starting from the end

print("Lastest released version\nGencode: ",GencodeVersion,"\nGenome: ",GenomeVersion,sep="")
RessourcesODIR="RESSOURCES/GENCODE/GencodeV"+GencodeVersion
GTF=RessourcesODIR+"/gencode.v"+GencodeVersion+".annotation.gtf"
Genome=RessourcesODIR+"/"+GenomeVersion+".genome.fa"


Hisat2IndexDir="RESSOURCES/GENCODE/GencodeV"+GencodeVersion+"/Hisat2IndexDir/"
Hisat2IndexRoot=Hisat2IndexDir+GenomeVersion+"_x_"+GencodeVersion
SSfile=Hisat2IndexDir+"SpliceSite.txt"
Hisat2ALN_DIR="ALN_FILES/HISAT2_"+GenomeVersion+"_x_"+GencodeVersion+"/"

FQ_DIR="/bioinfo/users/ltaing/ProteoGenomics/data/2000553/"

SAMPLES=["MB01","MB02","MB25"]

rule all:
	input: "blastp.outfmt6"

rule GetGTF:
	output: GTF
	shell:
		"mkdir --parents {RessourcesODIR} &&\
 wget \"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GencodeVersion}/gencode.v{GencodeVersion}.annotation.gtf.gz\"\
 --output-document {GTF}.gz &&\
 gunzip {GTF}.gz"

rule GetGenome:
	output: Genome
	shell: "mkdir --parents {RessourcesODIR} &&\
 wget \"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GencodeVersion}/{GenomeVersion}.genome.fa.gz\"\
 --output-document {Genome}.gz &&\
 gunzip {Genome}.gz"

rule BuildGenomeIndex:
	input: Genome
	output: Genome+".fai"
	conda: "envs/samtools.yaml"
	threads: 1
	shell: "samtools faidx -Sb {Genome}"

rule BuildHST2SpliceSite:
	input: GTF
	output: SSfile
	conda: "envs/Hisat2.yaml"
	shell: "mkdir --parents {Hisat2IndexDir} &&\
 hisat2_extract_splice_sites.py -v {GTF} > {SSfile}"

rule BuildHST2_INDEX:
	input: SSfile, Genome
	output: Hisat2IndexRoot+".1.ht2"
	conda: "envs/Hisat2.yaml"
	threads: 10
	shell: "hisat2-build -p {threads} \
 --ss {SSfile} {Genome} {Hisat2IndexRoot}"

rule Hisat2Map:
	input: index=Hisat2IndexRoot+".1.ht2", splicesite=SSfile,
	 R1=FQ_DIR+"{Sample}_R1.fastq.gz", R2=FQ_DIR+"{Sample}_R2.fastq.gz"
	output: SAM=Hisat2ALN_DIR+"{Sample}.sam",NSS=Hisat2ALN_DIR+"{Sample}.NovelSpliceSite.txt"
	conda: "envs/Hisat2.yaml"
	threads: 10
	shell: "mkdir --parents {Hisat2ALN_DIR};\
 hisat2 -x {Hisat2IndexRoot} -q --time --threads 10 \
 -1 {input.R1} -2 {input.R2} \
 --known-splicesite-infile {input.splicesite} \
 --novel-splicesite-outfile  {output.NSS} \
 --dta-cufflinks \
 -S {output.SAM}"

rule SamToBam:
	input: SAM=Hisat2ALN_DIR+"{Sample}.sam"
	output: BAM=Hisat2ALN_DIR+"{Sample}.bam"
	conda: "envs/samtools.yaml"
	threads: 1
	shell: "samtools view -Sb {input.SAM} > {output.BAM}"

rule SamToSortedBam:
	input: BAM=Hisat2ALN_DIR+"{Sample}.bam"
	output: SBAM=Hisat2ALN_DIR+"{Sample}.sorted.bam"
	conda: "envs/samtools.yaml"
	threads: 10
	shell: "samtools sort --threads {threads} -T {output.SBAM}.temp -o {output.SBAM} {input.BAM}"

rule FilterBamAndIndex:
	input: SBAM=Hisat2ALN_DIR+"{Sample}.sorted.bam"
	output: FSBAM=Hisat2ALN_DIR+"{Sample}.filtered.sorted.bam"
	shell: "samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 512 -F 2048 -q 60 {input.SBAM} > {output.FSBAM};\
 samtools index {output.FSBAM}"

rule IsoformsWithStrawberry:
	input: FSBAM=Hisat2ALN_DIR+"{Sample}.filtered.sorted.bam"
	output: "STRAWBERRY/{Sample}_strawberry_with_Ref.gtf"
	threads: 1
	shell: "mkdir --parents STRAWBERRY/{wildcards.Sample}/ ;\
 {Strawberry} --verbose --min-mapping-qual 60 -g {GTF} \
 --output-dir STRAWBERRY/{wildcards.Sample}/ {input.FSBAM} ;\
 mv STRAWBERRY/{wildcards.Sample}/assembled_transcripts.gtf {output};\
 rm --recursive --force STRAWBERRY/{wildcards.Sample}/"

rule MergeGTF:
	input: expand("STRAWBERRY/{Sample}_strawberry_with_Ref.gtf",Sample=SAMPLES)
	output: "STRAWBERRY/merge.ChrsOnly.gtf"
	conda: "envs/stringtie.yaml"
	shell: "stringtie \
 --merge \
 -o STRAWBERRY/merge.gtf \
 -c 2.5 \
 -m 300 \
 -T 1 \
 -f .001 \
 -i \
 {input};\
 grep -e '^\#' STRAWBERRY/merge.gtf > STRAWBERRY/merge.ChrsOnly.gtf;\
 grep -e '^chr' STRAWBERRY/merge.gtf >> STRAWBERRY/merge.ChrsOnly.gtf"

rule GTFToCDNA:
	input: GTF="STRAWBERRY/merge.ChrsOnly.gtf",GENOME=Genome,INDEX=Genome+".fai"
	output: "STRAWBERRY/merge.ChrsOnly.Cdna.fasta"
	conda: "envs/gffread.yaml"
	shell: "gffread {input.GTF} -w {output} -g {input.GENOME} "

rule LongOrfs:
	input: "STRAWBERRY/merge.ChrsOnly.Cdna.fasta"
	output: "STRAWBERRY/transdecoder/longest_orfs.pep"
	shell: "TransDecoder.LongOrfs -t STRAWBERRY/merge.ChrsOnly.fasta;\
 mkdir --parents STRAWBERRY/transdecoder/;\
 mv merge.ChrsOnly.fasta.transdecoder_dir/* STRAWBERRY/transdecoder/.;\
 rm --recursive --force merge.ChrsOnly.fasta.transdecoder_dir;\
 rm --recursive --force merge.ChrsOnly.fasta.transdecoder_dir.__checkpoints_longorfs;"

rule makeblastdb:
	input: "RESSOURCES/UNIPROT/UP000005640_9606.fasta"
	output: ["RESSOURCES/UNIPROT/UP000005640_9606.fasta.pin","RESSOURCES/UNIPROT/UP000005640_9606.fasta.phr","RESSOURCES/UNIPROT/UP000005640_9606.fasta.psq"]
	shell: "makeblastdb -in {input} -dbtype prot"

rule blastp:
	input: pep="STRAWBERRY/transdecoder/longest_orfs.pep", blastdb=["RESSOURCES/UNIPROT/UP000005640_9606.fasta.pin","RESSOURCES/UNIPROT/UP000005640_9606.fasta.phr","RESSOURCES/UNIPROT/UP000005640_9606.fasta.psq"]
	output: "blastp.outfmt6"
	shell: "blastp \
 -num_threads 2 \
 -query {input.pep} \
 -db RESSOURCES/UNIPROT/UP000005640_9606.fasta \
 -max_target_seqs 1 \
 -outfmt 6 \
 -evalue 1e-5 \
 > {output}"




rule MakeTemplate:
	output:
		"mqpar_1.6.2.6.Template_US.xml"
	shell:
		"mono {MQ} --create mqpar_1.6.2.6.Template_US.xml"

rule MakeScratchyMQPAR_XML:
	input: 
		template="mqpar_1.6.2.6.Template_US.xml"
	output:
		ScratchMQPAR_XML=ScratchMQPAR_Name
	run:
		UTF8_XML_PARSER=etree.XMLParser(encoding='UTF-8',remove_blank_text=True)
		tree=etree.parse(input.template,UTF8_XML_PARSER)
		
		Fasta = tree.xpath("/MaxQuantParams/fastaFiles/FastaFileInfo/fastaFilePath")[0]
		Fasta.text=config["MaxQuant"]["FASTA_DATABASE"]
		
		Session = tree.xpath("/MaxQuantParams/name")[0]
		Session.text="Job"+os.environ["PBS_JOBID"]
		
		Mods= tree.xpath("/MaxQuantParams/parameterGroups/parameterGroup/labelMods")[0]
		ZeMods=etree.Element("string")
		ZeMods.text="Arg10;Lys8"
		Mods.append(ZeMods)
		MaxLab= tree.xpath("/MaxQuantParams/parameterGroups/parameterGroup/maxLabeledAa")[0]
		MaxLab.text="3"
		Multiplicity= tree.xpath("/MaxQuantParams/parameterGroups/parameterGroup/multiplicity")[0]
		Multiplicity.text="2"
		fixedModifications=tree.xpath("/MaxQuantParams/parameterGroups/parameterGroup/fixedModifications/string")[0]
		fixedModifications.getparent().remove(fixedModifications)
		variableModifications=tree.xpath("/MaxQuantParams/parameterGroups/parameterGroup/variableModifications")[0]
		VarMod=etree.Element("string")
		VarMod.text="Carbamidomethyl (C)"
		variableModifications.append(VarMod)
		
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
					RmIndex="rm --recursive --verbose --force "+RmIndex
					os.system(RmIndex)
					RmFolder=line
					RmFolder=RmFolder.replace(".raw","")
					RmFolder="rm --recursive --verbose --force "+RmFolder
					os.system(RmFolder)
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
		"MaxQuantRapport.txt"
	threads: 24
	shell:
		"mono {MQ} {ScratchMQPAR_Name} > MaxQuantRapport.txt"
