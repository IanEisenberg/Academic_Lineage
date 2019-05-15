#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 22:39:25 2019

@author: ian
"""

import glob
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from os import path


def get_abstract_file(node, source='Data/*'):
    abstract = glob.glob(path.join(source, '%07d.abs' % node))[0]
    return abstract

def read_hepth_abstract(abstract, key='Abstract'):
    """ read HEP-Th abstract from SNAP
        source: https://snap.stanford.edu/data/cit-HepTh.html
    """
    out = {}
    f = open(abstract,'r').read().split('\\')
    # extra title and authors
    meta = f[2].split('\n')
    for line in meta:
        try:
            k,*v = line.split(': ')
            out[k.lower()] = ': '.join(v)
        except ValueError:
            pass
    # extract author
    abstract = f[4]
    abstract = abstract.replace('\n', ' ')
    out['abstract'] = abstract
    if key:
        out = out[key]
    return out
        
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
        self.graph = nx.from_edgelist(edge_list, create_using=nx.DiGraph)
        
    def get_lineage(self, node):
        """ 
        Outputs references and citations 
        for a node (paper) in the graph 
        """
        out_edges = self.graph.edges(nbunch=node)
        in_edges =  self.graph.in_edges(nbunch=node)
        return {'references': sorted([i[1] for i in out_edges]),
                'citations': sorted([i[0] for i in in_edges])}
    
    def get_lineage_abstract_files(self, node, source='Data/*'):
        """ 
        Gets abstract files for node
        """
        lineage = self.get_lineage(node)
        if type(source) == str:
            citation_abstracts = [get_abstract_file(i) for i in lineage['citations']]
            reference_abstracts = [get_abstract_file(i) for i in lineage['references']]
        return {'citations': citation_abstracts, 
                'references': reference_abstracts}
    
    def get_lineage_abstracts(self, node, source='Data/*'):
        abstracts = self.get_lineage_abstract_files(node, source)
        abstracts = {k: list(map(read_hepth_abstract, v)) for k,v in abstracts.items()}
        return abstracts


if __name__ == "__main__":
    edge_list = np.loadtxt('Data/cit-HepTh.txt.gz', dtype=int)
    graph = AncestryGraph()
    graph.load_from_edgelist(edge_list)
    node = 9611127
    lineage = graph.get_lineage(node)
    abstracts = graph.get_lineage_abstracts(node)
    f,axes = plt.subplots(2,2, figsize=(12,12))
    for i, key in enumerate(abstracts.keys()):
        # plot wordclouds
        ax = axes[i][0]
        abstracts_wordcloud(abstracts[key], ax=ax)
        ax.set_title(key.title(), fontsize=20, fontweight='bold')
        ax.axis('off')
        
        # plot citation graph
        # get subgraph of citations/references
        lineage_nodes = lineage[key]
        # show community one out
        additional_neighborhood = list(np.unique([i[1] for i in list(graph.graph.out_edges(nbunch=lineage_nodes))]))
        subgraph = graph.graph.subgraph(nodes=lineage_nodes+additional_neighborhood+[node])
        
        
        
        graph_pos = nx.spring_layout(subgraph, k=0.15)
        nx.draw_networkx_nodes(subgraph, pos=graph_pos,
                               nodelist=additional_neighborhood,
                               node_color='grey', 
                               alpha=.5,
                               node_size=100,
                               ax=axes[i][1])
        nx.draw_networkx_nodes(subgraph, pos=graph_pos,
                               nodelist=lineage_nodes,
                               node_color='b', node_size=250,
                               ax=axes[i][1])
        nx.draw_networkx_nodes(subgraph, pos=graph_pos,
                               nodelist=[node],
                               node_color='red', 
                               node_size=500,
                               ax=axes[i][1])
        if key=='references':
            nx.draw_networkx_edges(subgraph, pos=graph_pos,
                                   edgelist=subgraph.edges(nbunch=node),
                                   ax=axes[i][1])
        else:
            nx.draw_networkx_edges(subgraph, pos=graph_pos,
                                   edgelist=subgraph.in_edges(nbunch=node),
                                   ax=axes[i][1])
        axes[i][1].axis('off')
        
