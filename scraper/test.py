import scholarly
from IPython import embed

import paper_extraction, pdf_extraction, keywords, file_utility

def main():
	search_query = "bert nlp"

	print("Fetching paper matching: %s" % search_query)

	results = scholarly.search_pubs_query(search_query)
	paper = next(results)

	paper_data = paper_extraction.extract(paper, verbose=True)
	file_utility.save_json(paper_data, PAPER_DATA_PATH)

	embed()

def test_pdf_extraction():
	pdf_extraction.extract("https://arxiv.org/pdf/1811.06965")
	pass

def test_keywords():
	paper_data = file_utility.load_json(PAPER_DATA_PATH)

	abstracts = []
	for paper in paper_data["citing_papers"]:
		if "abstract" in paper:
			abstracts.append(paper["abstract"])

	keywords.extract_uncommon_words(abstracts[0])
	pass

PAPER_DATA_PATH = "./tmp/paper_data.json"

# main()
# test_pdf_extraction()
test_keywords()