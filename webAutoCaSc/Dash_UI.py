import copy
import io
import os
import tempfile

from flask import Flask, send_file

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from urllib.parse import unquote, quote
import pandas as pd
from AutoCaSc_core.AutoCaSc import AutoCaSc, AUTOCASC_VERSION
from dash_extensions import Download
from dash_extensions.snippets import send_data_frame

server = Flask(__name__)
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP],
                server=server)
app.title = "AutoCaSc"

url_bar_and_content_div = html.Div([
    dcc.Location(id='url', refresh=False),
    dcc.Store(id="query_memory"),
    dcc.Store(id="variant_queue"),
    dcc.Store(id="variant_memory"),
    dcc.Store(id="results_memory"),
    Download(id="download"),
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
                        dbc.Col(dbc.NavLink("Home", href='/'),
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
    dbc.Container(
        [
            html.Br(),
            html.H2("About"),
            html.P("AutoCaSc is a tool for quantifying the plausibility of candidate variants for Neurodevelopmental Disorders (NDD). AutoCaSc is intended to be used on impactful rare variants in a research setting. In its current version, 12 parameters are counted in, achieving a maximum of 15 points. User inputs are the identified variant in a standard HGVS/VCF format together with segregation aspects (de novo, recessive, dominant and X-chromosomal). We use the Ensembl REST API to annotate variant attributes (e.g. variant consequence, allele frequency from gnomAD, in silico predictions) and gene based scores dependent on inheritance mode (e.g. high Z-score is of relevant for de novo missense) from dbNSFP. Other attributes were previously labor intensive and predisposed to variability. These included important categories like expression in the nervous system, neuronal functions, co-expression and protein interactions, search for relevant literature, model organisms and observations in screening studies. As an objective approach we now searched a defined set of databases (GTEx, STRING, MGI, PubTator, PsyMuKB, DisGeNET) and generated empirical cut-offs for each category by comparing the respective readout between a manually curated list of known NDD genes from the SysID database and a list of genes not involved in NDDs."),
            html.Br(),
            html.P("Feel free to contact johann.lieberwirth@medizin.uni-leipzig.de in case you have further questions or in case you have found a bug.")
        ]
    )
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
    dbc.Button(id="download_button"),
    html.Div([
            navbar,
            dbc.Container([
                dbc.Card([
                    dbc.CardHeader(
                        dbc.Tabs(
                            id="card_tabs",
                            card=True,
                        )
                    ),
                    dbc.CardBody(
                        html.P(id="card_content", className="card_text")
                    )
                ])
            ])
        ]),
    html.Div(id="loading_output")
])


########## FRONTEND ##########
# basic string_formatting and initiating AutoCaSc intsances in order to check if input is ok
def parse_input_to_instances(input):
    variants = [input.split(",")[i].strip() for i in range(len(input.split(",")))]
    instances = [AutoCaSc(_variant, mode="web") for _variant in variants]
    return instances


def get_display_variant(_variant, n_chars=25):
    if len(_variant) < (n_chars + 5):
        return _variant
    else:
        return _variant[:n_chars] + "..."


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
            html.Br(),
            dbc.Container([
                dbc.Card([
                    dbc.CardHeader(
                        dbc.Row(
                            [
                                dbc.Tabs(
                                    tab_list,
                                    id="card_tabs",
                                    card=True,
                                    active_tab=initial_tab,
                                    # style={
                                    #     "margin-bottom": "0px !important"
                                    # }
                                ),

                            ],
                            className="align-items-baseline",
                            style={
                                "padding-bottom": "0 !important",
                                "margin-bottom": "0 !important",
                                "padding-left": "10px",
                                "padding-right": "10px",
                            },
                            justify="between",
                        ),
                    ),
                    dbc.CardBody(
                        html.P(id="card_content", className="card_text")
                    )
                ])
            ])
        ])
        return results_page


def input_ok(instances):
    for instance in instances:
        if instance.variant_format == "incorrect":
            return False
    return True


def get_badge(status_code):
    if status_code == 200:
        return None
    if status_code in [498, 496, 498]:
        return dbc.Badge(status_code, color="danger", className="mr-1")
    else:
        return dbc.Badge(status_code, color="warning", className="mr-1")



@app.callback(
        Output(f"navbar-collapse", "is_open"),
        [Input(f"navbar-toggler", "n_clicks")],
        [State(f"navbar-collapse", "is_open")])
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


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
        if "/search" in pathname and results_memory is None:
            print("search page")
            return search_page
        else:
            return landing_page
    else:
        print("results page")
        return get_results_page(results_memory)


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


@app.callback(
    Output("url", "pathname"),
    [Input("search_button", "n_clicks")],
    [State("variant_queue", "data"),
     State("inheritance_input", "value")]
)
def search_button_click(n_clicks, variant_queue, inheritance):
    if not n_clicks is None and not variant_queue is None and not inheritance is None:
        variants = [variant_queue.get("instances").get(_key).get("variant") for _key in variant_queue.get("instances").keys()]
        url_suffix = quote(f"/search/inheritance={inheritance}/variants={'%'.join(variants)}")
        return url_suffix
    else:
        raise PreventUpdate


@app.callback(
    Output("query_memory", "data"),
    Input("url", "pathname")
)
def interpret_url(pathname):
    if "search" in pathname:
        inheritance = unquote(pathname).split("inheritance=")[-1].split("/")[0]
        query_data = {}
        query_data["inheritance"] = inheritance
        return query_data
    raise PreventUpdate


@app.callback(
    Output("card_content", "children"),
    [Input("card_tabs", "active_tab")],
    [State("results_memory", "data")],
    State("card_content", "children")
)
def get_tab_card(active_tab, results_memory, old_card):
    cell_style = {
        "padding": "5px",
        "padding-left": "12px"
    }
    if active_tab == "overview_tab":
        overview_table_header = html.Thead(
            html.Tr(
                [
                    html.Th("Variant"),
                    html.Th("Candidate Score"),
                    html.Th("HGVS")
                ],
            )
        )
        overview_table_rows = []
        tooltips = []
        for i, _variant in enumerate(results_memory.get("instances").keys()):
            _instance_attributes = results_memory.get('instances').get(_variant)
            if _instance_attributes.get("gene_symbol"):
                _hgvs = _instance_attributes.get("gene_symbol")
                if _instance_attributes.get("hgvsc_change") is not None:
                    _hgvs += ":" + _instance_attributes.get("hgvsc_change")
                if _instance_attributes.get("hgvsp_change") is not None:
                    if ":" not in _hgvs:
                        _hgvs += ":"
                    _hgvs += _instance_attributes.get("hgvsp_change")
            else:
                _hgvs = None

            overview_table_rows.append(
                html.Tr(
                    [
                        html.Th([_variant, " ",
                                 get_badge(_instance_attributes.get("status_code"))],
                                style=cell_style),
                        html.Th(_instance_attributes.get("candidate_score_v1"),
                                style=cell_style),
                        html.Th(html.P(_hgvs,
                                       id=f"tooltip_target_{i}"),
                                style=cell_style)
                    ]
                )
            )
            tooltips.append(dbc.Tooltip(f"Transcript: {_instance_attributes.get('transcript')}",
                            target=f"tooltip_target_{i}"))
        overview_table = dbc.Table(
            [
                overview_table_header,
                html.Tbody(
                    overview_table_rows
                ),
            ],
            responsive=True,
            hover=True,
            striped=True,
        )
        tab_card_content = [
            html.Div(tooltips),
            dbc.Row(
                [
                    dbc.Col(html.H3("Overview on variants")),
                    dbc.Col(dbc.Button("Download",
                                       id="download_button",
                                       style={
                                           "margin-bottom": "10px",
                                           "margin-top": "0"
                                       }),
                            width="auto"
                            )
                ],
                justify="between",
            ),
            html.Br(),
            overview_table
        ]
        return tab_card_content
    else:
        tab_num = int(active_tab.split("_")[-1])
        _variant = list(results_memory.get("instances").keys())[tab_num]
        _instance_attributes = results_memory.get("instances").get(_variant)
        status_code = _instance_attributes.get("status_code")

        if len(results_memory.get("instances")) == 1:
            card_header = dbc.Row(
                    [
                        dbc.Col(html.H3(f"Variant: {get_display_variant(_variant)}"),
                                className="col-12 col-md-6"),
                        dbc.Col(
                            dbc.Row(
                                [
                                    dbc.Col(html.H3(f"Candidate Score: {_instance_attributes.get('candidate_score_v1')}")),
                                    dbc.Col(dbc.Button("Download",
                                                       id="download_button",
                                                       style={
                                                           "margin-bottom": "10px",
                                                           "margin-top": "0"
                                                       }),
                                            width="auto"),
                                ],
                            ),
                            className="col-12 col-md-6"
                        )
                    ]
                )
        else:
            card_header = dbc.Row(
                [
                    dbc.Col(html.H3(f"Variant: {get_display_variant(_variant)}"),
                            className="col-12 col-md-6"),
                    dbc.Col(html.H3(f"Candidate Score: {_instance_attributes.get('candidate_score_v1')}"),
                            className="col-12 col-md-6")
                ]
            )
        if status_code == 200:
            scoring_results = {
                "inheritance_score": "Inheritance",
                "gene_attribute_score": "Gene Attributes",
                "variant_score": "Variant Attributes",
                "literature_score": "Literature Plausibility",
            }

            parameters = {
                "impact": "Impact",
                "cadd_phred": "CADD phred",
                "oe_lof_interval": "o/e LoF",
                "oe_mis_interval": "o/e mis",
                "ac_hom": "gnomAD homozygous count",  # todo show gnomad hemi only when X
                "n_hemi": "gnomAD hemizygous count",
                "gerp_rs_rankscore": "GERP++ RS"
            }

            casc_table_header = html.Thead(
                                html.Tr(
                                    [
                                        html.Th("Subscore"),
                                        html.Th("")
                                    ]
                                )
                            )
            casc_table_rows = []
            for _subscore in scoring_results.keys():
                casc_table_rows.append(
                    html.Tr(
                        [
                            html.Th(scoring_results.get(_subscore),
                                    id=f"{_subscore}_description",
                                    scope="row",
                                    style=cell_style),
                            html.Td(_instance_attributes.get(_subscore),
                                    id=f"{_subscore}_explanation",
                                    style=cell_style)
                        ],
                    )
                )

            explanation_tooltips = [
                dbc.Tooltip(f"{_instance_attributes.get('explanation_dict').get('pli_z')}",
                            target="gene_attribute_score_explanation"),
                dbc.Tooltip(f"impact: {_instance_attributes.get('explanation_dict').get('impact')}, "
                            f"insilico: {_instance_attributes.get('explanation_dict').get('in_silico')}, "
                            f"conservation: {_instance_attributes.get('explanation_dict').get('conservation')}, "
                            f"frequency: {_instance_attributes.get('explanation_dict').get('frequency')}",
                            target="variant_score_explanation"),
                dbc.Tooltip(f"{_instance_attributes.get('explanation_dict').get('inheritance')}",
                            target="inheritance_score_explanation")
                # Todo literature score explanation
            ]
            description_tooltips = [
                dbc.Tooltip("Points attributed for inheritance & segregation",
                            target="inheritance_score_description"),
                dbc.Tooltip("Points attributed for gene constraint parameters.",
                            target="gene_attribute_score_description"),
                dbc.Tooltip("Points attributed for insilico predictions, conservation and allele frequency",
                            target="variant_score_description"),
                dbc.Tooltip("Points b for data in literature databases, animal models, expression pattern, "
                            "interaction networks",
                            target="literature_score_description")
            ]

            parameter_table_header = html.Thead(
                html.Tr(
                    [
                        html.Th("Parameter"),
                        html.Th("")
                    ],
                )
            )
            parameter_table_rows = []
            for _parameter in parameters.keys():
                parameter_table_rows.append(
                    html.Tr(
                        [
                            html.Th(parameters.get(_parameter),
                                    scope="row",
                                    style=cell_style),
                            html.Td(_instance_attributes.get(_parameter),
                                    style=cell_style)
                        ],
                    )
                )

            tab_card_content = [
                html.Div(description_tooltips + explanation_tooltips),
                card_header,
                html.Hr(),
                dbc.Row(
                    [
                        dbc.Col(dcc.Markdown(f"**Gene symbol:** {_instance_attributes.get('gene_symbol')}"),
                                className="col-12 col-md-6"),
                        dbc.Col(dcc.Markdown(f"**Transcript:** {_instance_attributes.get('transcript')}"),
                                className="col-12 col-md-6")
                    ],
                ),
                dbc.Row(
                    [
                        dbc.Col(dcc.Markdown(f"**HGVS:** {_instance_attributes.get('hgvsc_change')} "
                                             f"{_instance_attributes.get('hgvsp_change') or ''}"),
                                className="col-12 col-md-6"),
                        dbc.Col(dcc.Markdown(f"**VCF:** {_instance_attributes.get('vcf_string')}"),
                                className="col-12 col-md-6")
                    ],
                ),
                html.Br(),
                dbc.Row(
                    [
                        dbc.Col(
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
                                ),
                            className="col-12 col-md-6"
                        ),
                        dbc.Col(
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
                            ),
                            className="col-12 col-md-6"
                        ),
                    ],
                )
            ]
            return tab_card_content
        elif status_code == 201:
            return dbc.Alert("Error: The reference sequence does not match GRCh37!", color="warning")
        elif status_code == 201:
            return dbc.Alert("Error: The reference sequence does not match GRCh37!", color="warning")
        elif status_code == 400:
            return dbc.Alert("Error: Could not process variant. Please try VCF annotation or "
                             "HGVS annotation using HGNC gene symbol!",
                             color="danger")
        elif status_code == 496:
            return dbc.Alert("Error: The alternative sequence matches the GRCh37 reference sequence!", color="danger")
        elif status_code == 498:
            return dbc.Alert("Error: You have entered an intergenic variant.", color="danger")
        else:
            return dbc.Alert(f"Some error occured. Code {status_code}", color="warning")


########## BACKEND ##########
def score_variants(instances, inheritance):
    if not inheritance == "comphet":
        for _instance in instances:
            _instance.inheritance = inheritance
            if _instance.__dict__.get("status_code") == 200:
                _instance.calculate_candidate_score()
    else:
        variant_gene_df = pd.DataFrame()
        for i, _instance in enumerate(instances):
            _instance.inheritance = inheritance
            if _instance.__dict__.get("status_code") == 200:
                variant_gene_df.loc[i, "variant"] = _instance.__dict__.get("variant")
                variant_gene_df.loc[i, "gene_symbol"] = _instance.__dict__.get("gene_symbol")
                variant_gene_df.loc[i, "instance"] = _instance

        instances = []
        for _gene in variant_gene_df.gene_symbol.unique():
            df_chunk = variant_gene_df.loc[variant_gene_df.gene_symbol == _gene].reset_index(drop=True)
            if len(df_chunk) != 2:
                raise IOError(f"combination of comphet variants unclear for {_gene}")
                #todo error output when corresponding comphet variant not found
            for i in range(2):
                _instance = copy.deepcopy(df_chunk.loc[abs(1-i), "instance"])
                _instance.other_autocasc_obj = df_chunk.loc[abs(0-i), "instance"]
                _instance.calculate_candidate_score()
                instances.append(_instance)
    return instances


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
    Output("results_memory", "data"),
    Input("variant_memory", "data"),
    Input("query_memory", "data")
)
def calculate_results(variant_memory, query_memory):
    inheritance = None
    if query_memory:
        inheritance = query_memory.get("inheritance")
    if variant_memory is not None and inheritance is not None:
        instances = dict_to_instances(variant_memory)
        instances = score_variants(instances, inheritance)
        return store_instances(instances)
    raise PreventUpdate


@app.callback(
    Output("download", "data"),
    Input("download_button", "n_clicks"),
    State("results_memory", "data")
)
def download_button_click(n_cklicks, results_memory):
    if not n_cklicks:
        raise PreventUpdate
    df = pd.DataFrame()
    for i, _variant in enumerate(results_memory.get("instances").keys()):
        _instance = results_memory.get("instances").get(_variant)
        df.loc[i, "variant"] = _instance.get("variant")
        df.loc[i, "gene_symbol"] = _instance.get("gene_symbol")
        df.loc[i, "transcript"] = _instance.get("transcript")
        df.loc[i, "hgvsc"] = _instance.get("hgvsc_change")
        df.loc[i, "hgvsp"] = _instance.get("hgvsp_change")
        df.loc[i, "impact"] = _instance.get("impact")
        df.loc[i, "inheritance"] = _instance.get("inheritance")
        df.loc[i, "candidate_score"] = _instance.get("candidate_score_v1")
        df.loc[i, "literature_plausibility"] = _instance.get("literature_score")
        df.loc[i, "inheritance_score"] = _instance.get("inheritance_score")
        df.loc[i, "variant_attribute_score"] = _instance.get("variant_score")
        df.loc[i, "gene_attribute_score"] = _instance.get("gene_attribute_score")

    data = io.StringIO()
    df.to_csv(data, sep="\t", decimal=",")
    data.seek(0)
    return dict(content=data.getvalue(), filename="AutoCaSc_results.tsv")


if __name__ == '__main__':
    app.run_server(debug=True,
                   dev_tools_hot_reload=True,
                   host='0.0.0.0',
                   port=5000)