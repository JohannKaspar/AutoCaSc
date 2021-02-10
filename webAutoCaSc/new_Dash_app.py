import copy
import time

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from urllib.parse import unquote, quote
import pandas as pd
from AutoCaSc_core.AutoCaSc import AutoCaSc, AUTOCASC_VERSION


app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

url_bar_and_content_div = html.Div([
    dcc.Location(id='url', refresh=False),
    dcc.Store(id="query_memory"),
    dcc.Store(id="variant_queue"),
    dcc.Store(id="variant_memory"),
    dcc.Store(id="results_memory"),
    dcc.Interval(id="check_for_data_interval", n_intervals=0, interval=200, disabled=True),
    html.Div(id='page-content',
             style={"width":"100%",
                    "color": "#000000"})
])


# navbar = dbc.NavbarSimple(
#     children=[
#         dbc.NavItem(dbc.NavLink("Page 1", href="/")),
#         dbc.NavItem(dbc.NavLink("About", href='/about')),
#         dbc.NavItem(dbc.NavLink("another link", href='/about')),
#     ],
#     brand="AutoCaSc",
#     brand_href="#",
#     color="dark",
#     dark=True,
# )

navbar = dbc.Navbar(
    dbc.Container(
        [
            html.A(
                dbc.NavbarBrand("AutoCaSc",
                                style={
                                    "font-size": "1.8em"
                                }),
                href="/",
            ),
            # dbc.Col(dbc.NavLink("About", href='/about'),
            #                 width="auto"),
            # dbc.Col(dbc.NavLink("another link", href='/about'),
            #         width="auto"),
            dbc.NavbarToggler(id="navbar-toggler"),
            dbc.Collapse(
                dbc.Row(
                    [
                        dbc.Col(dbc.NavLink("About", href='/about'),
                                width="auto"),
                        dbc.Col(dbc.NavLink("another link", href='/about'),
                                width="auto"),
                        # dbc.Col(
                        #     dbc.Input(type="search", placeholder="Search variant")
                        # ),
                        # dbc.Col(
                        #     dbc.Button(
                        #         "Search", color="primary", className="ml-2"
                        #     ),
                        #     # set width of button column to auto to allow
                        #     # search box to take up remaining space.
                        #     width="auto",
                        # ),
                    ],
                    no_gutters=True,
                    className="ml-auto flex-nowrap mt-3 mt-md-0",
                    align="center",
                ),
                id="navbar-collapse",
                navbar=True,
            )
        ],
    fluid=True
    ),
    color="dark",
    dark=True,
)

# navbar = dbc.Navbar(
#     [
#         dbc.Col(dbc.NavbarBrand("AutoCaSc", href="/"),
#                 sm=2,
#                 md=1),
#         dbc.Col(dbc.NavLink("About", href='/about'),
#                 width="auto"),
#         dbc.Col(dbc.NavLink("another link", href='/about'),
#                 width="auto"),
#         dbc.Col(dbc.Row(
#             [
#                 dbc.Col(dbc.Input(type="search", placeholder="Search variant")),
#                 dbc.Col(dbc.Button("Search", color="primary", className="ml-2 mr-2"),
#                         width="auto")
#             ],
#             no_gutters=True,
#             className="ml-auto flex-nowrap mt-3 mt-md-0",
#             align="center",
#             ))
#     ],
#     color="dark",
#     dark=True,
#     #align="center",
# )

variant_input_card = dbc.FormGroup(
    [
        dbc.Label("Enter a variant here"),
        dbc.Input(
            type="text",
            id="variant_input",
            placeholder="e.g. X:12345:T:C",
            # value="1:7725246:G:A",  #todo delete this
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
        dbc.Spinner(fullscreen=True)
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
    html.Div(style={"height": "15vh"}),

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
                ),
                html.Br(),
                dbc.Row(
                    [
                        dbc.Col(
                            dbc.Button(
                                "Start search",
                                color="primary",
                                id="search_button"
                            )
                        ),
                        dbc.Col(
                            html.Div(
                                id="spinner_card"
                            )
                        )
                    ]
                )
            ],
        ),

])



# index layout
app.layout = url_bar_and_content_div

# "complete" layout
app.validation_layout = html.Div([
    url_bar_and_content_div,
    about_page,
    landing_page,
    search_page,
    results_page_clear,
    html.Div(id="loading_output")
])


# basic string_formatting and initiating AutoCaSc intsances in order to check if input is ok
def parse_input_to_instances(input):
    variants = [input.split(",")[i].strip() for i in range(len(input.split(",")))]
    instances = [AutoCaSc(_variant, mode="web") for _variant in variants]
    return instances


def get_display_variant(_variant):
    if len(_variant) < 25:
        return _variant
    else:
        return _variant[:17] + "..."

def input_ok(instances):
    for instance in instances:
        if instance.variant_format == "incorrect":
            return False
    return True


def get_parameter(parameter, dict):
    return dict.get(parameter) or "-"

def get_results_page(results_memory): #todo add error handling
    if not results_memory:
        return results_page_clear
    else:
        tab_list = []
        initial_tab = "tab_0"
        tab_num = 0
        if len(results_memory.get("instances").items()) > 1:
            tab_list.append(dbc.Tab(
                label="Overview",
                tab_id="overview_tab"
            ))
            initial_tab = "overview_tab"
        for _variant, _instance_attributes in results_memory.get("instances").items():
            tab_list.append(dbc.Tab(label=get_display_variant(_variant), tab_id=f"tab_{tab_num}"))
            tab_num += 1

        results_page = html.Div([
            navbar,
            dbc.Container([
                dbc.Card([
                    dbc.CardHeader(
                        dbc.Tabs(
                            tab_list,
                            id="card_tabs",
                            card=True,
                            active_tab=initial_tab
                        )
                    ),
                    dbc.CardBody(
                        html.P(id="card_content", className="card_text")
                    )
                ])
            ])
        ])
        return results_page

@app.callback(
    Output("card_content", "children"),
    [Input("card_tabs", "active_tab")],
    [State("results_memory", "data")]
)
def get_tab_card(active_tab, results_memory):
    #todo add gene name etc
    if active_tab == "overview_tab":
        variants = list(results_memory.get("instances").keys())
        candidate_scores = [results_memory.get("instances").get(_variant).get("candidate_score_v1")
                            for _variant in variants]
        tab_card_content = [
            html.H3("Overview on variants"),
            dbc.Table.from_dataframe(pd.DataFrame(
                {
                    "Variant": variants,
                    "Candidate Score": candidate_scores
                }),
                striped=True,
                bordered=True,
                hover=True,
                responsive=True
            )
        ]
    else:
        tab_num = int(active_tab.split("_")[-1])
        _variant = list(results_memory.get("instances").keys())[tab_num]
        _instance_attributes = results_memory.get("instances").get(_variant)

        parameters = {
            "cadd_phred": "CADD phred",
            "oe_lof_interval": "o/e LoF",
            "oe_mis_interval": "o/e mis",
            "ac_hom": "gnomAD homozygous count",  # todo show gnomad hemi only when X
            "n_hemi": "gnomAD hemizygous count",
            "gerp_rs_rankscore": "GERP++ RS"
        }
        scoring_results = {
            "inheritance_score": "Inheritance",
            "gene_attribute_score": "Gene Attributes",
            "variant_score": "Variant Attributes",
            "literature_score": "Literature Plausibility",
            #"candidate_score_v1": "Sum"
        }

        cell_style = {
            "padding": "5px",
            "padding-left": "12px"
        }

        casc_table_header = html.Thead(
                            html.Tr(
                                [
                                    html.Th("Subscore",
                                        className="col-6"),
                                    html.Th("",
                                        className="col-6")
                                ],
                                className="d-flex"
                            )
                        )
        casc_table_rows = []
        for _subscore in scoring_results.keys():
            casc_table_rows.append(
                html.Tr(
                    [
                        html.Th(scoring_results.get(_subscore),
                                scope="row",
                                style=cell_style,
                                className="col-6"),
                        html.Td(_instance_attributes.get(_subscore),
                                style=cell_style,
                                className="col-6")
                    ],
                    className="d-flex"
                )
            )

        parameter_table_header = html.Thead(
            html.Tr(
                [
                    html.Th("Parameter",
                            className="col-6"),
                    html.Th("",
                            className="col-6")
                ],
                className="d-flex"
            )
        )
        parameter_table_rows = []
        for _parameter in parameters.keys():
            parameter_table_rows.append(
                html.Tr(
                    [
                        html.Th(parameters.get(_parameter),
                                scope="row",
                                style=cell_style,
                                className="col-6"),
                        html.Td(_instance_attributes.get(_parameter),
                                style=cell_style,
                                className="col-6")
                    ],
                    className="d-flex"
                )
            )

        tab_card_content = [
            dbc.Row(
                [
                    dbc.Col(html.H3(f"Variant: {get_display_variant(_variant)}")),
                    dbc.Col(html.H3(f"Candidate Score: {_instance_attributes.get('candidate_score_v1')}"))
                ]
            ),
            html.Br(),
            dbc.Row(
                [
                    dbc.Col(dcc.Markdown(f"**Gene symbol:** {_instance_attributes.get('gene_symbol')}\n")),
                    dbc.Col(dcc.Markdown(f"**Transcript:** {_instance_attributes.get('transcript')}"))
                ],
            ),
            # html.Br(),
            dbc.Row(
                dbc.Table(
                    [
                        casc_table_header,
                        html.Tbody(
                            casc_table_rows
                        ),
                    ],
                    responsive=True,
                    hover=True,
                    striped=True,
                    #style={"width": "100%"}
                )
            ),
            dbc.Row(
                dbc.Table(
                    [
                        parameter_table_header,
                        html.Tbody(
                            parameter_table_rows
                        ),
                    ],
                    responsive=True,
                    hover=True,
                    striped=True,
                    style={"width": "100%"}
                )
            )
        ]
    return tab_card_content


def score_variants(instances, inheritance):
    if not inheritance == "comphet":
        for _instance in instances:
            _instance.inheritance = inheritance
            _instance.calculate_candidate_score()
    else:
        variant_gene_df = pd.DataFrame()
        for i, _instance in enumerate(instances):
            _instance.inheritance = inheritance
            variant_gene_df.loc[i, "variant"] = _instance.__dict__.get("variant")
            variant_gene_df.loc[i, "gene_symbol"] = _instance.__dict__.get("gene_symbol")
            variant_gene_df.loc[i, "instance"] = _instance

        instances = []
        for _gene in variant_gene_df.gene_symbol.unique():
            df_chunk = variant_gene_df.loc[variant_gene_df.gene_symbol == _gene].reset_index(drop=True)
            if len(df_chunk) != 2:
                raise IOError(f"combination of comphet variants unclear for {_gene}")
            for i in range(2):
                _instance = copy.deepcopy(df_chunk.loc[abs(1-i), "instance"])
                _instance.other_autocasc_obj = df_chunk.loc[abs(0-i), "instance"]
                _instance.calculate_candidate_score()
                instances.append(_instance)

    return instances

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
    [Input("variant_memory", "data")],
    [State("query_memory", "data")]
)
def get_results(variant_memory, query_memory):
    if variant_memory is not None and query_memory is not None:
        inheritance = query_memory.get("inheritance")
        instances = dict_to_instances(variant_memory)
        instances = score_variants(instances, inheritance)
        return store_instances(instances)
    raise PreventUpdate


# "#spinner
# @app.callback(,
#               [Input('url', 'pathname'),
#                Input("results_memory", "data")])
# def trigger_spinner(pathname, memory_data):
#     ctx = dash.callback_context
#     if "pathname" in ctx.triggered[0]['prop_id'] and "search" in pathname:
#         time.sleep(30)
#         return ""
#     return """


# Index callbacks
@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname'),
               Input("results_memory", "data")])
def display_page(pathname, results_memory):
    #todo add possibility to directly enter variants in url
    ctx = dash.callback_context
    if "pathname" in ctx.triggered[0]['prop_id']:
        if pathname == "/about":
            return about_page
        if pathname == "/results":
            if results_memory is not None: #return empty site just in case
                raise PreventUpdate
        if "/search" in pathname:
            return search_page
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
    [Input("variant_input", "n_blur")],
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

def instances_to_dict(instance):
    instance_dict = {}
    for key, value in instance.__dict__.items():
        if any([x in type(value).__name__ for x in ["int", "float", "bool", "NoneType", "str", "dict", "list"]]):
            instance_dict[key] = value
        else:
            instance_dict[key] = instances_to_dict(value)
    return instance_dict

def store_instances(instance_list, code_key="variant"):
    instance_dicts = [instances_to_dict(_instance) for _instance in instance_list]
    return {"instances":{_instance_dict.get(code_key): _instance_dict for _instance_dict in instance_dicts}}


@app.callback(
    Output("variant_memory", "data"),
    [Input("variant_queue", "data")],
    [State("variant_memory", "data")]
)
def retrieve_variant_data(variant_queue, variant_memory):
    if variant_queue:
        if variant_memory is not None:
            if variant_queue.get("instances").keys() == variant_memory.get("instances").keys():
                return store_instances(dict_to_instances(variant_memory))
        instances = dict_to_instances(variant_queue)
        for _instance in instances:
            if not _instance.data_retrieved:
                _instance.retrieve_data()
        # [_instance.retrieve_data() for _instance in instances if not _instance.data_retrieved]
        return store_instances(instances)
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
    if not n_clicks is None and not variant_queue is None and not inheritance is None:
        variants = [variant_queue.get("instances").get(_key).get("variant") for _key in variant_queue.get("instances").keys()]
        url_suffix = quote(f"/search/inheritance={inheritance}/variants={'%'.join(variants)}")

        query_data = {}
        query_data["inheritance"] = inheritance
        spinner = dbc.Spinner(
                                html.Div(id="loading_output"),
                                fullscreen=True,
                                spinner_style={
                                    "visibility": "hidden"
                                }
                            )

        return url_suffix, False, query_data
    else:
        raise PreventUpdate


if __name__ == '__main__':
    app.run_server(debug=False,
                   dev_tools_hot_reload=False)