import dash
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from dash.dependencies import Input, Output

import numpy as np
import pandas as pd
# custom
from utils import get_abstract_file, read_hepth_abstract

# style dictionaries
table_cell_style = {"fontFamily": "Arial",
                    "size": 10, 
                     'textAlign': 'left',
                     'whiteSpace': 'no-wrap',
                     'overflow': 'hidden',
                     'textOverflow': 'ellipsis'}

# other setup
test_df = pd.DataFrame(np.random.rand(10,4))
colors = {
    'background': '#111111',
    'text': '#7FDBFF'
}

app = dash.Dash()

app.layout = html.Div(
    [
    # Title - Row
    html.Div(
        [
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
                'A tool for exploring citation networks around scientific papers',
                style={'font-family': 'Helvetica',
                       "font-size": "120%",
                       "width": "80%"},
                className='eight columns',
            ),
        ],
        className='row'
    ),

    # selectors
    html.Div([
        html.Label('Paper Selection'),
        dcc.Dropdown(
            className='four columns',
            id='paper-select',
            options=[
                {'label': 'Paper 1', 'value': 9611127},
                {'label': 'Paper 2', 'value': 9611126},
            ],
            value=9611127,
        )],
        className='row'
    ),

    # abstracts
    html.Div(
        [
            html.Div(
                className="four columns",
                style={'padding': 15},
                children = [
                # select paper
                html.H5(
                    id='title-div',
                    style={'textAlign': 'center'}),
                html.H6(
                    id='authors-div',
                    style={'textAlign': 'center' }),
                html.Div(
                    id='abstract-div',
                    className='collapsible'
                )
            ]),
            html.Div([
                dash_table.DataTable(
                    id='reference-abstracts',
                    columns=[{"name": i, "id": i, 'deletable': True} for i in test_df], 
                    row_selectable="multi",
                    selected_rows=[0],
                    style_cell = table_cell_style,
                    data=test_df.to_dict("rows")
                    )
                ], className="four columns"),
            html.Div([
                dash_table.DataTable(
                    id='citation-abstracts',
                    columns=[{"name": i, "id": i, 'deletable': True} for i in test_df], 
                    row_selectable="multi",
                    selected_rows=[0],
                    style_cell = table_cell_style,
                    data=test_df.to_dict("rows")
                    )
                ], className="four columns")

        ],
        style={'height': "33vh"},
        className='row'
    ),
    # in-focus abstract and word clouds
        html.Div(
        [
            html.Div(
                className="four columns",
                children = [
                # select paper
                html.H4(
                    id='title-div2',
                    style={'textAlign': 'center'}),
                html.H6(
                    id='authors-div2',
                    style={'textAlign': 'center' }),
                html.Div(
                    id='abstract-div2',
                    className='collapsible')
            ]),
            html.Div(
                className="four columns",
                children=html.Div([
                    dcc.Graph(
                        id='reference-word-cloud',
                    )
                ])
            ),
            html.Div(
                className="four columns",
                children=html.Div([
                    dcc.Graph(
                        id='citation-word-cloud',
                    ),
                ])
            )
        ],
        style={'height': "33vh"},
        className='row'
    )
])

# set up call backs
@app.callback(
    [Output('abstract-div', 'children'),
     Output('title-div', 'children'),
     Output('authors-div', 'children')],
    [Input('paper-select', 'value')]
)
def update_primary_abstract(input_value):
    abstract_data = read_hepth_abstract(get_abstract_file(input_value), None)
    return ("Abstract: " + abstract_data['Abstract'], 
            abstract_data['Title'], 
            abstract_data['Authors'])

@app.callback(
    [Output('abstract-div2', 'children'),
     Output('title-div2', 'children'),
     Output('authors-div2', 'children')],
    [Input('paper-select', 'value')]
)
def update_target_abstract(input_value):
    abstract_data = read_hepth_abstract(get_abstract_file(input_value), None)
    return abstract_data['Abstract'], abstract_data['Title'], abstract_data['Authors']



if __name__ == '__main__':
    app.run_server(port=8051, debug=True)

