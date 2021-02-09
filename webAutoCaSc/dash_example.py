import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from urllib.parse import unquote, quote
import pandas as pd

DBC_DOCS = "https://dash-bootstrap-components.opensource.faculty.ai/"
DBC_GITHUB = "https://github.com/facultyai/dash-bootstrap-components"

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

url_bar_and_content_div = html.Div([
    dcc.Location(id='url', refresh=False),
    dcc.Store(id="query_memory"),
    dcc.Store(id="variant_queue"),
    dcc.Store(id="variant_memory"),
    dcc.Store(id="results_memory"),
    dcc.Interval(id="check_for_data_interval", n_intervals=0, interval=200, disabled=True),
    html.Div(id='page-content')
])


navbar = dbc.NavbarSimple(
    children=[
        dbc.NavItem(dbc.NavLink("GitHub", href=DBC_GITHUB)),
        dbc.DropdownMenu(
            nav=True,
            in_navbar=True,
            label="Menu",
            children=[
                dbc.DropdownMenuItem("Entry 1", href="https://google.com"),
                dbc.DropdownMenuItem("Entry 2", href="/test"),
                dbc.DropdownMenuItem(divider=True),
                dbc.DropdownMenuItem("A heading", header=True),
                dbc.DropdownMenuItem(
                    "Entry 3", href="/external-relative", external_link=True
                ),
                dbc.DropdownMenuItem("Entry 4 - does nothing"),
            ],
        ),
    ],
    brand="Dash Bootstrap Components",
    brand_href=DBC_DOCS,
    sticky="top",
)

variant_input_card = dbc.FormGroup(
    [
        dbc.Label("Enter a variant here"),
        dbc.Input(
            type="text",
            id="variant_input",
            placeholder="e.g. X:12345:T:C",
            value="1:7725246:G:A",  #todo delete this
            autoFocus=True
        ),
        dbc.FormText("Although HGVS works as well, we recommend using VCF format."),
    ]
)

misc_input_card = dbc.FormGroup(
    [
        dbc.RadioItems(
            id="inheritance_input",
            options=[
                {"label": "de novo", "value": "de_novo"},
                {"label": "inherited dominant", "value": "ad_inherited"},
                {"label": "homozygous recessive", "value": "homo"},
                {"label": "X-linked", "value": "x_linked"},
                {"label": "compound heterozygous", "value": "comphet"},
            ],
            inline=True
        ),
        # html.Br(),
        # dbc.RadioItems(
        #     id="assembly_input",
        #     options=[
        #         {"label": "GRCh37 (hg19)", "value": "GRCh37"},
        #         {"label": "GRCh38 (hg38)", "value": "GRCh38"},
        #     ],
        #     inline=True,
        #     value="GRCh37",
        # )
    ]
)

search_page = html.Div([
    navbar,
    dbc.Container([
        dbc.Spinner() #todo format spinner
        # html.P(
        #     ""
        # )
    ])
])

results_page_clear = html.Div([
    navbar,
    dbc.Container([
        html.H3(
            "X:12345:C:T"
        ),
        html.H2(
            "Candidate Score: 9"
        ),
        html.P(
            "Gene: XYZ"
        ),
        html.P(
            "HGVSC: ENST100001023:c.213C>T"
        ),
        html.P(
            "HGVSP: XYZ:p.Thr70Asn"
        )
        # html.P(
        #     ""
        # )
    ],
        #fluid=True
    )
],

)

about_page = html.Div([
    navbar,
    html.H2('Impressum'),
    html.P("This page belongs to Johann Lieberwirth. Blablabla®")
])


landing_page = html.Div([
    navbar,
    dbc.Container(
        [
            dbc.Row(
                dbc.Col(
                    [
                        html.H2('Enter your variant below'),
                        variant_input_card,
                        misc_input_card
                    ]
                ),
                align="baseline",
                justify="center"
            ),
            html.Br(),
            dbc.Row(
                dbc.Col(
                    dbc.Button(
                        "Start search",
                        color="primary",
                        id="search_button"
                    )
                )
            )
        ]
    )
])

table = html.Div(
    [
        html.H2("Table"),
        dbc.Table(
            [
                html.Thead(
                    html.Tr(
                        [
                            html.Th("#"),
                            html.Th("First name"),
                            html.Th("Last name"),
                        ]
                    )
                ),
                html.Tbody(
                    [
                        html.Tr(
                            [
                                html.Th("1", scope="row"),
                                html.Td("Tom"),
                                html.Td("Cruise"),
                            ]
                        ),
                        html.Tr(
                            [
                                html.Th("2", scope="row"),
                                html.Td("Jodie"),
                                html.Td("Foster"),
                            ]
                        ),
                        html.Tr(
                            [
                                html.Th("3", scope="row"),
                                html.Td("Chadwick"),
                                html.Td("Boseman"),
                            ]
                        ),
                    ]
                ),
            ],
            responsive=True,
            striped=True,
            hover=True,
        ),
    ]
)


app.layout = html.Div(
    [
        navbar,
        dbc.Container(
            [
                html.Br(),
                table,
                html.Div(style={"height": "200px"}),
            ]
        ),
    ]
)




@app.callback(
    Output("collapse", "is_open"),
    [Input("collapse-button", "n_clicks")],
    [State("collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

def score_variants(instances, inheritance):
    if not inheritance == "comphet":
        for _instance in instances:
            _instance.inheritance = inheritance
            _instance.calculate_candidate_score()
        return instances
    else:  #todo add beaviour for comphets
        pass



# @app.callback(
#     Output("results_memory", "data"),
#     Input("check_for_data_interval", "n_intervals"),
#     [State("variant_memory", "data"),
#      State("query_memory", "data")]
# )
# def get_results(_, variant_memory, query_memory):
#     if variant_memory is None:
#         raise PreventUpdate
#     else:
#         instances = dict_to_instances(variant_memory)
#         instances = score_variants(instances)
#         return {"result": "nice"}

@app.callback(
    Output("results_memory", "data"),
    [Input("variant_memory", "data"),
     Input('url', 'pathname')],
    [State("query_memory", "data")]
)
def get_results(variant_memory, pathname, query_memory):
    if "/search/" in pathname:
        if variant_memory is not None and query_memory is not None:
            inheritance = query_memory.get("inheritance")
            instances = dict_to_instances(variant_memory)
            instances = score_variants(instances, inheritance)
            return instances_to_dict(instances)
    raise PreventUpdate

# Index callbacks
@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname'),
               Input("results_memory", "data")])
def display_page(pathname, results_memory):
    ctx = dash.callback_context
    if "pathname" in ctx.triggered[0]['prop_id'] and "search" not in pathname:
        if pathname == "/about":
            return about_page, False
        if pathname == "/results":
            # return get_results_page(None)  #return empty site just in case
            return get_results_page(None)  #return empty site just in case
        # if "/search/" in pathname:
        #     return get_results_page(pathname, variant_memory, inheritance="dnovo")
        else:
            return landing_page
    else:
        return get_results_page(results_memory)


@app.callback(
        Output(f"navbar-collapse", "is_open"),
        [Input(f"navbar-toggler", "n_clicks")],
        [State(f"navbar-collapse", "is_open")])
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("variant_input", "valid"),
    Output("variant_input", "invalid"),
    Output("variant_queue", "data"),
    [Input("variant_input", "n_blur")],  #Todo better use "debounce" instead of n_blur?
    [State("variant_input", "value")]
)
def check_user_input(_, user_input):
    if not user_input == None:
        variant_instances = parse_input_to_instances(user_input)
        if input_ok(variant_instances):
            variant_queue = {}
            variant_queue["instances"] = {_instance.__dict__.get("variant"): _instance.__dict__ for _instance in variant_instances}
            return True, False, variant_queue
        else:
            return False, True, None
    raise PreventUpdate


def dict_to_instances(dict):
    instances = []
    try:
        for _variant in dict.get("instances").keys():
            _instance = dict.get("instances").get(_variant)
            instance = AutoCaSc(_variant, mode="web")
            if _instance.get("data_retrieved"):
                for _key in _instance.keys():
                    instance.__dict__[_key] = _instance.get(_key)
            instances.append(instance)
        return instances
    except AttributeError:
        return None

def instances_to_dict(instance_list, key="variant"):
    return {"instances":{_instance.__dict__.get(key): _instance.__dict__ for _instance in instance_list}}


@app.callback(
    Output("variant_memory", "data"),
    [Input("variant_queue", "data")],
)
def retrieve_variant_data(variant_queue):
    if variant_queue:
        instances = dict_to_instances(variant_queue)
        [_instance.retrieve_data() for _instance in instances]
        return instances_to_dict(instances)
    else:
        raise PreventUpdate


@app.callback(
    [Output("url", "pathname"),
     Output("check_for_data_interval", "disabled"),
     Output("query_memory", "data")],
    [Input("search_button", "n_clicks")],
    [State("variant_queue", "data"),
     State("inheritance_input", "value")]
)
def search_button_click(n_clicks, variant_queue, inheritance):
    #todo add spinner in case of waiting
    if not n_clicks is None and not variant_queue is None and not inheritance is None:
        variants = [variant_queue.get("instances").get(_key).get("variant") for _key in variant_queue.get("instances").keys()]
        url_suffix = quote(f"/search/inheritance={inheritance}/variants={'%'.join(variants)}")

        query_data = {}
        query_data["inheritance"] = inheritance
        return url_suffix, False, query_data
    else:
        raise PreventUpdate



if __name__ == "__main__":
    app.run_server(port=8888, debug=True)