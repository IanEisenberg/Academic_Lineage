import urllib
import subprocess

def extract(url):
	download_pdf(url)
	extract_text()

def extract_text():
	command = "pdftotext %s %s" % (TEMP_PAPER_PDF_PATH, TEMP_PAPER_TEXT_PATH)
	subprocess.call([command], shell=True)

def download_pdf(url):
	pdf_data = urllib.request.urlopen(url).read()
	file = open(TEMP_PAPER_PDF_PATH, "wb")
	file.write(pdf_data)
	file.close()

TEMP_PAPER_PDF_PATH = "./tmp/temp-paper.pdf"
TEMP_PAPER_TEXT_PATH = "./tmp/temp-paper.txt"