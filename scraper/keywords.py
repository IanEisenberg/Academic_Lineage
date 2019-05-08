import file_utility
import nltk

from IPython import embed


COMMON_WORD_PATH = "./config/english-words-top-1000.txt"
common_words = set(file_utility.load_text(COMMON_WORD_PATH).split("\n"))

def extract_uncommon_words(text):
	tokens = set(x.lower() for x in nltk.word_tokenize(text))
	embed()
