#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 22:39:25 2019

@author: ian
"""

import networkx
import numpy as np

class AncestryGraph():
    def __init__(self, edge_list):
        self.graph = networkx.from_edgelist(edge_list, 
                                            create_using=networkx.DiGraph)
    
    def get_lineage(self, node):
        out_edges = self.graph.edges(nbunch=node)
        in_edges =  self.graph.in_edges(nbunch=node)
        return {'ancestry': [i[1] for i in out_edges],
                'children': [i[0] for i in in_edges]}
        

if __name__ == "__main__":
    edge_list = np.loadtxt('Data/cit-HepTh.txt.gz')
    graph = AncestryGraph(edge_list)
    print(graph.get_lineage(9301122.0))
