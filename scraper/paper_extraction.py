import scholarly
from IPython import embed
from colorama import Fore, Back, Style
from bs4 import BeautifulSoup

import urllib, re

def extract(paper, fetch_abstract=True, verbose=False):
	cited_count = paper.citedby

	citing_paper_gen = paper.get_citedby()

	citing_papers = []

	index = 0
	while True:
		if verbose:
			print("Fetching citing paper %i of %i" % (index, cited_count))

		try:
			citing_paper = next(citing_paper_gen)
			citing_papers.append(paper_output(citing_paper, fetch_abstract=fetch_abstract))

		except StopIteration:
			break

		index += 1

	output = paper_output(paper)
	if index == 0:
		print(output)

	output["citing_papers"] = citing_papers
	return output


def paper_output(paper, fetch_abstract=False):
	output = {
		"title": paper.bib.get("title"),
		"author": paper.bib.get("author"),
		"meta_url": paper.bib.get("url"),
		"paper_url": paper.bib.get("eprint"),
	}

	try:
		output["cited_count"] = paper.citedby
	except AttributeError:
		output["cited_count"] = 0

	if fetch_abstract and output.get("meta_url"):
		extracted = extract_from_website(output["meta_url"])

		if "abstract" in extracted:
			output["abstract"] = extracted["abstract"]

	return output

def extract_from_website(url):
	extracted = {}

	if re.fullmatch(ARXIV_PATTERN, url):
		html = urllib.request.urlopen(url).read()
		soup = BeautifulSoup(html, "html.parser")

		abstracts = soup.select(".abstract")
		if len(abstracts):
			extracted["abstract"] = abstracts[0].text

	return extracted

ARXIV_PATTERN = ".+arxiv\.org\/abs\/.+"
