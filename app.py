#!/usr/bin/env python3


# dynamic controls: https://community.plot.ly/t/dynamic-controls-and-dynamic-output-components/5519
# chains callbacks: https://community.plot.ly/t/order-of-chained-callbacks-when-app-initialises-or-call-only-some-callbacks-on-start-up/6015
#!/usr/bin/env python
import biopython_parser as bp
import dash
import dash_table
import dash_core_components as dcc
import dash_html_components as html
from plotly.tools import mpl_to_plotly
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# custom
from utils import (get_abstract_file, 
                   read_hepth_abstract, 
                   AncestryGraph)
from plot_utils import abstracts_wordcloud, fig_to_uri

# style dictionaries
table_cell_style = {"fontFamily": "Arial",
                    "size": 10, 
                     'textAlign': 'left',
                     'whiteSpace': 'no-wrap',
                     'overflow': 'hidden',
                     'textOverflow': 'ellipsis'}

# set up temp data

# components
# Title - Row
banner = html.Div([
            html.H1(
                'Academic Lineage',
                style={'font-family': 'Helvetica',
                       "margin-top": "25",
                       "margin-bottom": "0"},
                className='eight columns',
            ),
            html.Img(
                src="https://cayley.io/img/carousel/cayley.png",
                className='two columns',
                style={
                    'height': '9%',
                    'width': '9%',
                    'float': 'right',
                    'position': 'relative',
                    'padding-top': 10,
                    'padding-right': 0
                },
            ),
            html.P(
                'A tool for exploring the citation context around scientific papers',
                style={'font-family': 'Helvetica',
                       "font-size": "120%",
                       "width": "80%"},
                className='eight columns',
            ),
        ], className='row banner')

selectors = html.Div([
                html.Div([
                        html.Label('Search Field'),
                        dcc.Dropdown(
                                id='field-select',
                                options=[
                                        {'label': 'Keywords', 'value': ''},
                                        {'label': 'Paper Title', 'value': 'title'},
                                        {'label': 'Paper Author', 'value': 'author'},
                                        {'label': 'PubMed ID', 'value': 'pmid'}
                                        ],
                                value='',
                                ),
                        ], className='one columns'),
            html.Div([
                        html.Label('0 Results', id='index-label'),
                        dcc.Dropdown(
                                id='index-select',
                                options=[
                                        {'label': '1', 'value': '0'},
                                        ],
                                value='0',
                                ),
                        ], className='one columns'),
            html.Div([
                html.Label('Search Query'),
                dcc.Input(
                    type='text',
                    id='paper-select',
                    value='brain PLOS', # 'Enter a paper title or PMID here',
                    className='twelve columns'
                ),
            ], className='two columns'),
            html.Div([
                html.Label('Lineage Selector'),
                dcc.RadioItems(
                    id='lineage-select',
                    options=[
                        {'label': 'Citations', 'value': 'citations'},
                        {'label': 'References', 'value': 'references'},
                    ],
                    value='references',
                    labelStyle={'display': 'inline-block'})
                ], className='two columns'),
            html.Div([
                html.Label('Number: 0', id='number-label'),
                ], className='two columns'),
        ], className='row selectors')

# paper, references and citations
main_row = html.Div( [
            html.Div(
                className="four columns subcontext",
                style={'padding': 15},
                children = [
                # select paper
                html.P(
                    'Hover to scroll and read abstract',
                    id='header-primary',
                    style={'textAlign': 'center', 'font-weight': "bold"}),
                html.Div(
                    id='abstract-primary',
                    className='collapsible'
                )
            ]),
            html.Div([
                html.H6(
                    id='lineage-header',
                    style={'textAlign': 'center' }),
                dcc.RadioItems(
                    id='lineage-list',
                    style={'max-height': '30vh',
                          'overflow': 'auto'}
                    ),
                ], className="four columns"),
            html.Div(
                className="four columns subcontext",
                children = [
                # select paper
                html.P(
                    'Hover to scroll and read abstract',
                    id='header-target',
                    style={'textAlign': 'center', 'font-weight': "bold"}),
                html.Div(
                    id='abstract-target',
                    className='collapsible'
                )
            ]),

        ],
        style={'height': "33vh"},
        className='row main')

# context row
context_row = html.Div([
                html.Div(
                    className="four columns",
                    children=html.Div([
                        html.Img(
                            id='word-cloud')
                        ])
                    ),
                html.Div(
                    className='eight columns context',
                    children=[
                            html.Div(
                                className="six columns subcontext",
                                children = [
                                html.H4('Context within paper',
                                    id='context-title',
                                    style={'textAlign': 'center' }),
                                html.Div(
                                    id='context-text',
                                    className='collapsible subcontext-pane')
                                ]
                            ),
                            html.Div(
                                className="six columns subcontext",
                                children = [
                                html.H6('Links',
                                    id='explore-title',
                                    style={'textAlign': 'center' }),
                                html.A(
                                    html.Button( "Link to Paper on GScholar"),
                                    href='https://google.com',
                                    target="_blank",
                                    id='link-primary'
                                    ),
                                html.A(
                                    html.Button("Link to Lineage Paper on GScholar"),
                                    href='https://google.com',
                                    target="_blank",
                                    id='link-target'
                                    )
                                ]),
                        ]
                    ),
                ],
                style={'height': "33vh"},
                className='row'
            )
                
                
# set up dashboard
app = dash.Dash()
app.layout = html.Div(
    [
    banner,
    selectors,
    main_row,
    context_row,
    # Hidden div inside the app that stores the intermediate value
    html.Div(id='data_store', style={'display': 'none'}, children = bp.get_dash_data(None, index=0, field='')[0])   
    ],
    style={'backgroundColor':'#fffcef'}
)

# set up call backs
# update data to be used by rest of the components
@app.callback([Output('data_store', 'children'),
               Output('index-select', 'options'),
               Output('index-label', 'children')],
              [Input('paper-select', 'n_submit'),
               Input('index-select', 'value')],
              [State('paper-select', 'value'),
              State('field-select', 'value')])
def update_dash_data(ns1, index, query, field):
    index = int(index)
    if ns1 == 0:
        data, ids = bp.get_dash_data(None, index=index, field=field)
    else:
        data, ids = bp.get_dash_data(query, index=index, field=field)
    
    options = [{'label': str(i+1), 'value': str(i)} for i in range(len(ids))]
    return (data, options, '%s Results' % str(len(options)))



@app.callback(
    [Output('abstract-primary', 'children'),
     Output('link-primary', 'href')],
    [Input('paper-select', 'n_submit'), 
    Input('data_store', 'children')]
)
def update_primary_abstract(id, stored_data):
    if type(stored_data) == str:
        abstract_data = bp.load_dash_json(stored_data)
        authors = ', '.join(abstract_data['authors'])
        if len(authors) > 150:
            authors = authors[:140]+'...'
            
        to_fill = dcc.Markdown('###### %s (%s) \n ** %s ** \n\n Abstract: %s' % (abstract_data['title'], 
                                                                        abstract_data['date'],
                                                                        authors,
                                                                        abstract_data['abstract']))
        url_fill = '+'.join(abstract_data['title'].split(' '))
        url = 'https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=' + url_fill + '&btnG='
        return (to_fill, url)
    

@app.callback(
     [Output('abstract-target', 'children'),
      Output('link-target', 'href')],
     [Input('lineage-list', 'value'),
      Input('data_store', 'children'),
      Input('lineage-select', 'value')]
 )
def update_target_abstract(PMID, stored_data, selector):
    stored_data = bp.load_dash_json(stored_data)
    if stored_data['PMID'] is not None:
        if PMID in stored_data['lineage'][selector].keys():
            abstract_data = stored_data['lineage'][selector][PMID]
            authors = ', '.join(abstract_data['authors'])
            if len(authors) > 150:
                authors = authors[:140]+'...'
            to_fill = dcc.Markdown('###### %s (%s) \n ** %s ** \n\n Abstract: %s' % (abstract_data['title'], 
                                                                        abstract_data['date'],
                                                                        authors,
                                                                        abstract_data['abstract']))
            url_fill = '+'.join(abstract_data['title'].split(' '))
            url = 'https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=' + url_fill + '&btnG='
            return (to_fill, url)
        else:
            return ('No related paper selected', 'https://scholar.google.com/')
    else:
            return ('No related paper selected', 'https://scholar.google.com/')

    
        
# update citaiton/reference lists
def get_list_element(x, i):
    text = '%s (%s)' % (x['title'], x['date'])
    return text

@app.callback(
    [Output('lineage-list', 'options'),
     Output('number-label', 'children')],
    [Input('paper-select', 'n_submit'), 
    Input('data_store', 'children'),
    Input('lineage-select', 'value')]
)
def update_network_lists(id, stored_data, selector):
    stored_data = bp.load_dash_json(stored_data)
    if stored_data['PMID'] is not None:
        to_list = list(stored_data['lineage'][selector].values())
        titles = [get_list_element(x, i) for i,x in enumerate(to_list)]
        pmids = [x['PMID'] for x in to_list]
        out = []
        for i, t in zip(pmids, titles):
            out.append({'label': t, 'value': i})
        return (out, 'Number of Papers: %s' % len(titles))


# # update wordclouds
@app.callback(
    [Output('word-cloud', 'src')],
    [Input('paper-select', 'n_submit'), 
    Input('data_store', 'children')]
)
def update_wordcloud(id, stored_data):
    stored_data = bp.load_dash_json(stored_data)
    if stored_data['PMID'] is not None:
        citation_abs = [i['abstract'] for i in stored_data['lineage']['citations'].values()]
        reference_abs = [i['abstract'] for i in stored_data['lineage']['references'].values()] 
        abstracts = citation_abs + reference_abs
        if abstracts:
            wordcloud = abstracts_wordcloud(abstracts)
            return [fig_to_uri(wordcloud, transparent=True)]
        else:
            return [None]
    else:
            return [None]

# update context
@app.callback(
     [Output('context-text', 'children')],
     [Input('lineage-list', 'value'),
      Input('data_store', 'children'),
      Input('lineage-select', 'value')]
 )
def update_context(PMID, stored_data, selector):
    stored_data = bp.load_dash_json(stored_data)
    paper = stored_data['paper']
    if paper is not None and selector=='references' and PMID in paper['refs_to_paragraphs'].keys():
        references = paper['references']
        reference_id = [i+1 for i,r in enumerate(references) if PMID == r['pmid_cited']]
        if len(reference_id)>=1:
            reference_id = reference_id[0]
        else:
            reference_id = -1
        paragraph_ids = paper['refs_to_paragraphs'][PMID]
        paragraph_texts = [paper['paragraph'][i]['text'] for i in paragraph_ids]
        paragraph_texts = ['(Par %s) %s' % (paragraph_ids[i]+1,s) for i,s in enumerate(paragraph_texts)]
        joined = '\n'.join(paragraph_texts)
        parts = joined.split('[%s]'%reference_id)
        to_return = [parts[0]]
        for i in range(1, len(parts)):
            to_return.append(html.Span('[%s]' % reference_id, style={'color':'teal'}))
            to_return.append(parts[i])
        return [to_return]
    else:
        return ['']

 
if __name__ == '__main__':
    app.run_server(port=8051, debug=True)

