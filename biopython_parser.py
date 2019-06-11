from Bio import Entrez
from collections import defaultdict
import pubmed_parser as pp
import json

# Lineage helpers
def get_empty_info():
    out = {'PMID': None,
           'authors': ['Authors: '],
           'abstract': '',
           'title': 'Title: ',
           'date': 'No Date',
           'paper': None
           }
    return out

def get_article_info(entry, get_paper=False):
    """ 
    Extracts article info into datastructure for further use 
    """
    out = get_empty_info()
    # get info from Medline
    out['PMID'] = str(entry['MedlineCitation']['PMID'])
    article = entry['MedlineCitation']['Article']
    out['title'] = article['ArticleTitle']
    # check for authors
    if 'AuthorList' in article.keys():
        authorlist = article['AuthorList']
        out['authors'] = ['%s %s' % (i.get('ForeName','First:NA'), i.get('LastName', 'Last: NA')) for i in authorlist]
    else:
        out['authors'] = "No Authors Listed"
    # check for abstract
    if 'Abstract' in article:
        out['abstract'] = article['Abstract']['AbstractText']
        if type(out['abstract']) == list:
            out['abstract'] = ''.join(out['abstract'])
    else:
        out['abstract'] = 'No Abstract Found'
    # get date
    if 'ArticleDate' in article.keys() and len(article['ArticleDate'])>0:
        out['date'] = article['ArticleDate'][0]['Year']
    else:
        out['date'] = 'No Date Found'
    out['keywords'] = article.get('KeywordList','')
    # get paper if possible
    if get_paper:
        PMCID = PMID_to_PMCID(out['PMID'])
        if PMCID is not None:
            paper_info = download_paper(PMCID)
            out['paper_info'] = paper_info
    return out

def get_articles_info(id_list, get_paper=False):
    """
    Extracts info from id list
    """
    f = lambda x: get_article_info(x, get_paper)
    links_info = []
    if len(id_list) > 0:
        result = fetch_details(id_list)
        links_info = list(map(f, result))
    return links_info
    
def get_lineage_info(lineage_ids, get_paper=False):
    f = lambda x: get_articles_info(x, get_paper)
    return {k: {kk:vv for kk,vv in zip(v,f(v))} for k,v in lineage_ids.items()}

# test for problems
#get_paper=True
#for k,v in lineage_ids.items():
#    result = fetch_details(v)
#    for entry in result:
#        f(entry)
        

def get_dash_data(search_query, index=0, **kwargs):
    if search_query is None or search_query == '':
        return json.dumps(get_empty_info()), []
    id_list= search(search_query, **kwargs)['IdList']
    if len(id_list)>0:
        first_id = id_list[min(index, len(id_list)-1)]
    else:
        return json.dumps({}), []
    abstract_info = get_articles_info([first_id])[0]       
    PMCID = PMID_to_PMCID(first_id)
    paper_info = download_paper(PMCID)
    lineage_ids = get_lineage_ids([first_id])
    lineage_info = get_lineage_info(lineage_ids, get_paper=False)
    out = abstract_info
    out['paper'] = paper_info
    out['lineage'] = lineage_info
    out['similar'] = get_similar_ids(first_id)
    return json.dumps(out), id_list

def load_dash_json(dash_json):
    return json.loads(dash_json)

# Biopython helpers
def get_pretty_article_info(id_list, db='pubmed'):
    """ 
    Returns nicely formatted title, abstract, etc.
    """
    ids = ','.join(id_list)
    handle = Entrez.efetch(db=db, id=ids, retmode='text', rettype='abstract')
    return handle.read()

def search(query, **kwargs):
    """ 
    Generic search term to get articles 
    Lookup: https://dataguide.nlm.nih.gov/edirect/esearch.html
    """
    kwargs['sort'] = kwargs.get('sort', 'relevance')
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed', 
                            retmax='100',
                            retmode='xml', 
                            term=query,
                            **kwargs)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list, db='pubmed', subset='PubmedArticle'):
    """ 
    Given a set of ids, return pubmed article info
    """
    ids = ','.join(id_list)
    Entrez.email = 'ianeisenberg90@gmail.com'
    handle = Entrez.efetch(db=db,
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results[subset]

def get_links_ids(id_list, db='pubmed', linkname=''):
    """
    Return articles related to input id list
    See https://www.ncbi.nlm.nih.gov/pmc/tools/cites-citedby/ for more info
    Args:
        id_list: list of pubmed ids
        linkname: determines type of link. Left empty, all will be returned. Options:
            "pubmed_pubmed": returns similar articles retrieved using a word weight algorithm
            "pubmed_pubmed_citedin": articles that cite the articles in id_list
            "pubmed_pubmed_refs": articles that the articles in id_list reference
        
    """
    ids = ','.join(id_list)
    links = Entrez.elink(db=db, id=ids, linkname=linkname)
    records = Entrez.read(links)
    to_return = {}
    for record in records[0]['LinkSetDb']:
        link_list = [i['Id'] for i in record['Link']]
        linkname = record['LinkName']
        to_return[linkname] = link_list
    return to_return

def get_similar_ids(id_list):
    link_ids = get_links_ids(id_list)
    similar_ids = link_ids['pubmed_pubmed']
    return similar_ids
    
def get_lineage_ids(id_list):
    link_ids = get_links_ids(id_list)
    citation_ids = []
    reference_ids = []
    if 'pubmed_pubmed_citedin' in link_ids.keys():
        citation_ids = link_ids['pubmed_pubmed_citedin']
    if 'pubmed_pubmed_refs' in link_ids.keys():
        reference_ids = link_ids['pubmed_pubmed_refs']
    return {'citations': citation_ids,
           'references': reference_ids}

# extract info from paper
def PMID_to_PMCID(PMID):
    """ Tries to find a PMCID from a PMID"""
    handle = Entrez.elink(dbfrom="pubmed", db="pmc", linkname="pubmed_pmc", id=PMID, retmode="text")
    result = Entrez.read(handle)[0]['LinkSetDb']
    if result:
        return result[0]['Link'][0]['Id']
    else:
        return None

def download_paper(PMCID, get_paragraphs=True):
    if PMCID is None:
        return None
    handle = Entrez.efetch(db="pmc", rettype="full", retmode="xml", id=PMCID)
    xml=handle.read()
    try:
        refs = pp.parse_pubmed_references(xml)
    except AttributeError: # failed to find paper:
        return None
    if refs is not None:
        refs_dict = {i['ref_id']: (i['pmid_cited'] if i['pmid_cited'] else None) for i in refs}
        paragraphs = pp.parse_pubmed_paragraph(xml)
        references = pp.parse_pubmed_references(xml)
        ref_to_paragraphs = defaultdict(list)
        for i,paragraph in enumerate(paragraphs):
            if 'reference_ids' in paragraph:
                paragraph['reference_ids'] = [refs_dict.get(i, None) for i in paragraph['reference_ids']]
                for ref in paragraph['reference_ids']:
                    if ref:
                        ref_to_paragraphs[ref].append(i)
        return {'paragraph': paragraphs,
                'references': references,
               'refs_to_paragraphs': ref_to_paragraphs}
    else:
        return None

if __name__ == "__main__":
    # PMID = 19777066
    id_list = search('Regulation of Adipose Tissue Stromal Cells Behaviors by Endogenic Oct4 Expression Control')['IdList']
    PMCID = PMID_to_PMCID(id_list[0])
    paper_info = download_paper(PMCID)
    lineage_ids = get_lineage_ids(id_list)
    lineage_info = get_lineage_info(lineage_ids, get_paper=True)
    dash_info = load_dash_json(get_dash_data(id_list[0]))