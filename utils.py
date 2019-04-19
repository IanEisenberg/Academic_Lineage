#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 22:39:25 2019

@author: ian
"""

import glob
import matplotlib.pyplot as plt
import networkx
import numpy as np
from os import path
from wordcloud import WordCloud

def read_hepth_abstract(abstract):
    """ read HEP-Th abstract from SNAP
        source: https://snap.stanford.edu/data/cit-HepTh.html
    """
    f = open(abstract,'r').read()
    abstract = f.split('\\')[4]
    abstract = abstract.replace('\n', ' ')
    return abstract

def abstracts_wordcloud(abstracts, ax=None):
    """ take a list of abstract texts and output a wordcloud """ 
    wordcloud = WordCloud(max_font_size=50, 
                           max_words=100, 
                           background_color="white").generate(' '.join(abstracts))
    if ax is None:
        plt.figure(figsize=(12,8));
        plt.imshow(wordcloud, interpolation="bilinear")
        plt.axis("off")
    else:
        ax.imshow(wordcloud, interpolation="bilinear")
        
class AncestryGraph():
    """ 
    Loads citation graphs from different sources 
    and provides utility functions 
    """
    def __init__(self):
        self.graph = None
    
    def load_from_edgelist(self, edge_list):
        # the edges are directed. Edge from node i -> j means i cited j
        # so the out_degree of each node is the number of references
        # the in_degree is the number of citations to it
        self.graph = networkx.from_edgelist(edge_list, create_using=networkx.DiGraph)
        
    def get_lineage(self, node):
        """ 
        Outputs references (ancestors) and citations (children) 
        for a node (paper) in the graph 
        """
        out_edges = self.graph.edges(nbunch=node)
        in_edges =  self.graph.in_edges(nbunch=node)
        return {'ancestors': sorted([i[1] for i in out_edges]),
                'children': sorted([i[0] for i in in_edges])}
    
    def get_lineage_abstract_files(self, node, source='Data/*'):
        """ 
        Gets abstract files for node
        """
        lineage = self.get_lineage(node)
        if type(source) == str:
            child_abstracts = [glob.glob(path.join(source, '%07d.abs' % i))[0] for i in lineage['children']]
            ancestor_abstracts = [glob.glob(path.join(source, '%07d.abs' % i))[0] for i in lineage['ancestors']]
        return {'children': child_abstracts, 
                'ancestors': ancestor_abstracts}
    
    def get_lineage_abstracts(self, node, source='Data/*'):
        abstracts = self.get_lineage_abstract_files(node, source)
        abstracts = {k: list(map(read_hepth_abstract, v)) for k,v in abstracts.items()}
        return abstracts
    
if __name__ == "__main__":
    edge_list = np.loadtxt('Data/cit-HepTh.txt.gz', dtype=int)
    graph = AncestryGraph()
    graph.load_from_edgelist(edge_list)
    node = 9611127
    abstracts = graph.get_lineage_abstracts(node)
    
    f,axes = plt.subplots(2,1, figsize=(12,12))
    for ax in axes:
        ax.axis('off')
    for key, ax in zip(abstracts.keys(), axes):
        abstracts_wordcloud(abstracts[key], ax=ax)
        ax.set_title(key.title(), fontsize=20, fontweight='bold')
