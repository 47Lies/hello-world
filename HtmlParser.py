from bs4 import BeautifulSoup
from urllib import request
import re
DD=BeautifulSoup(request.urlopen('https://www.gencodegenes.org/releases/current.html'))
version_HTML=str(DD.find("h1"))
#'<h1>Release 28 (GRCh38.p12)</h1>' something like that
pre.GencodeVersion=re.sub("^[^ ]+ ","",version_HTML) #discard anything before the first space included starting from the begining
GencodeVersion=re.sub(" .*$","",pre.GencodeVersion) #discard anything after the first space included starting from the end


pre.GenomeVersion=re.sub("^[^(]+\(","",version_HTML) #discard anything before the first parenthesis starting from the begining
re.sub("\).*$","",pre.GenomeVersion) # discard anything after the the first parenthesis starting from the end
