import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html

import flask
import pandas as pd
import time
import os

server = flask.Flask('app')
server.secret_key = os.environ.get('secret_key', 'secret')

app = dash.Dash('app', server=server)

app.scripts.config.serve_locally = False

app.layout = html.Div([
    html.H1('This is a test'),
    html.A(
        html.Img(src="https://picsum.photos/200/300",
                 # "https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png",
                 height="30px"
                 ),
        href="https://creativecommons.org/licenses/by-nc-sa/4.0/",
        target="_blank"),
], className="container")

if __name__ == '__main__':
    app.run_server()