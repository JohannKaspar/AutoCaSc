import copy
import io
import os
import tempfile
from statistics import mean
import time

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
app.title = "webAutoCaSc"

url_bar_and_content_div = html.Div([
    dcc.Location(id='url', refresh=False),
    dcc.Store(id="query_memory"),
    dcc.Store(id="variant_queue_input"),
    dcc.Store(id="variant_queue_url"),
    dcc.Store(id="variant_memory"),
    dcc.Store(id="results_memory"),
    Download(id="download"),
    html.Div(id='page-content',
             style={"width":"100%",
                    "color": "#000000",
                    })
])

navbar = html.Div(
    [
        dbc.Navbar(
            dbc.Container(
                [
                    html.A(
                        dbc.NavbarBrand("webAutoCaSc",
                                        style={
                                            "font-size": "1.8em"
                                        }),
                        href="/",
                    ),
                    dbc.NavbarToggler(id="navbar-toggler"),
                    dbc.Collapse(
                        dbc.Row(
                            [
                                dbc.Col(dbc.NavLink("About", href='/about', style={"color": "#ffffff"}),
                                        width="auto"),
                                dbc.Col(dbc.NavLink("FAQ", href='/faq', style={"color": "#ffffff"}),
                                        width="auto"),
                                dbc.Col(dbc.NavLink("Impressum", href='/impressum', style={"color": "#ffffff"}),
                                        width="auto"),
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
            fixed="top",
            #sticky="top"
        ),
    html.Div(style={"height": "80px"})
    ]
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

footer = html.Div(
    [
        dbc.Navbar(
            dbc.Container(
                [
                    html.A(html.Img(src="https://mirrors.creativecommons.org/presskit/buttons/88x31/svg/by-nc-sa.eu.svg",
                                    height="30px"),
                           href="https://creativecommons.org/licenses/by-nc-sa/4.0/",
                           target="_blank"),
                    dbc.NavLink("Github", href='https://github.com', target="_blank", style={"color": "#ffffff"}),
                    dbc.NavLink("Human Genetics Leipzig",
                                href='https://www.uniklinikum-leipzig.de/einrichtungen/humangenetik',
                                target="_blank",
                                style={"color": "#ffffff"}),
                    dbc.NavLink("Our Manuscript",
                                href="https://www.biorxiv.org",
                                target="_blank",
                                style={"color": "#ffffff"})
                ],
            ),
            color="dark",
            dark=True,
            fixed="bottom"
        ),
        html.Div(style={"height": "70px"})
    ]
)

variant_input_card = dbc.FormGroup(
    [
        dbc.Input(
            type="text",
            id="variant_input",
            placeholder="e.g. X:12345:T:C",
            autoFocus=True
        ),
        #dbc.FormText(""),
        dcc.Markdown('''Although HGVS works as well, we recommend using VCF format.\t
                     Examples: [11:94730916:A:C](/search/inheritance%3Dde_novo/variants%3D11%3A94730916%3AA%3AC), 
                     [X:101409056:A:C](/search/inheritance%3Dx_linked/variants%3DX%3A101409056%3AA%3AC), 
                     [ENST00000378402.5:c.4966G>A](/search/inheritance%3Dhomo/variants%3DENST00000378402.5%3Ac.4966G%3EA)
                     ''',
                     style={"font-size": "12px",
                            "margin-top": "10px"})
    ]
)

misc_input_card = dbc.FormGroup(
    [
        dbc.RadioItems(
            id="inheritance_input",
            options=[
                {"label": "De novo", "value": "de_novo"},
                {"label": "Inherited dominant", "value": "ad_inherited"},
                {"label": "Homozygous recessive", "value": "homo"},
                {"label": "X-linked", "value": "x_linked"},
                {"label": "Compound heterozygous", "value": "comphet"},
                {"label": "Unknown", "value": "unknown"}
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

landing_page = html.Div([
    navbar,
    html.Div(style={"height": "10vh"}),
    dbc.Container(
            [
                dbc.Row(
                    dbc.Col(
                        [
                            dcc.Markdown("""# Welcome to **webAutoCaSc**,\n
#### a webinterface for the automatic CaSc classification of research candidate variants in neurodevelopmental disorders."""),
                            # html.Hr(),
                            dcc.Markdown("Enter your variant of interest and presumed inheritance mode here:"),
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
                        )
                    ]
                )
            ],
        # style={"min-height": "calc(89vh - 135px)"}
        ),
    # html.Div(style={"height": "1vh"}),
    footer
])


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
        style={"height": "90vh"}
    ),
    footer
],

)

citations = dcc.Markdown("""
---
1. [__VEP__](https://grch37.ensembl.org/info/docs/tools/vep/index.html): McLaren, W. et al. The Ensembl Variant Effect Predictor. Genome Biol 17, 122 (2016).\n
2. [__gnomAD__](https://gnomad.broadinstitute.org): Karczewski, K. J. et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 581, 434–443 (2020).\n
3. [__GTEx__](http://www.gtexportal.org/home/index.html): Consortium, T. Gte. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369, 1318–1330 (2020).\n
4. [__STRING__](https://string-db.org): Szklarczyk, D. et al. STRING v11: protein-protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets. Nucleic Acids Res 47, D607–D613 (2019).\n
5. [__MGI__](http://www.informatics.jax.org): Bult, C. J. et al. Mouse Genome Database (MGD) 2019. Nucleic Acids Res 47, D801–D806 (2019).\n
6. [__PubTator__](https://www.ncbi.nlm.nih.gov/research/pubtator/): Wei, C.-H., Allot, A., Leaman, R. & Lu, Z. PubTator central: automated concept annotation for biomedical full text articles. Nucleic Acids Res 47, W587–W593 (2019).\n
7. [__PsyMuKB__](http://www.psymukb.net): Lin, G. N. et al. PsyMuKB: An Integrative De Novo Variant Knowledge Base for Developmental Disorders. Genomics Proteomics Bioinformatics 17, 453–464 (2019).\n
8. [__DisGeNET__](https://www.disgenet.org): Piñero, J. et al. DisGeNET: a comprehensive platform integrating information on human disease-associated genes and variants. Nucleic Acids Res 45, D833–D839 (2017).\n
---
""")

about_page = html.Div([
    navbar,
    dbc.Container(
        [
            dbc.Row([
                html.H2("About"),
                dbc.Button("DE", id="about_language_button")
            ],
            justify="between",
            style={
                "padding-left": "18px",
                "padding-right": "18px",
            }),
            html.Br(),
            html.Div([
                html.P("AutoCaSc is a tool for quantifying the plausibility of candidate variants for Neurodevelopmental Disorders (NDD). AutoCaSc is intended to be used on impactful rare variants in a research setting. In its current version, 12 parameters are counted in, achieving a maximum of 15 points. User inputs are the identified variant in a standard HGVS/VCF format together with segregation aspects (de novo, recessive, dominant and X-chromosomal). We use the Ensembl VEP REST API (1) to annotate variant attributes (e.g. variant consequence, allele frequency from gnomAD, in silico predictions) and gene based scores dependent on inheritance mode (e.g. high Z-score is of relevant for de novo missense) from dbNSFP. Other attributes were previously labor intensive and predisposed to variability. These included important categories like expression in the nervous system, neuronal functions, co-expression and protein interactions, search for relevant literature, model organisms and observations in screening studies. As an objective approach we now searched a defined set of databases (gnomAD (2), GTEx (3), STRING (4), MGI (5), PubTator (6), PsyMuKB (7), DisGeNET (8)) and generated empirical cut-offs for each category by comparing the respective readout between a manually curated list of known NDD genes from the SysID database (9) and a list of genes not involved in NDD.",
                       style={"text-align": "justify"}),
                html.Br(),
                html.P("Feel free to contact johann.lieberwirth@medizin.uni-leipzig.de or rami.aboujamra@medizin.uni-leipzig.de in case you have further questions or in case you have found a bug.",
                       style={"text-align": "justify"})
            ],
            id="about_text"),
            html.Br(),
            citations
        ],
        # style={"min-height": "calc(99vh - 150px)"}
    ),
    # html.Div(style={"height": "1vh"}),
    footer
    ])

faq_page = html.Div([
    navbar,
    dbc.Container(
        [
            dcc.Markdown("""
            FAQ
            ---
            __What is AutoCaSc?__  
            The AutoCaSc tool systematically evaluates the plausibility of variants in genes not yet associated with human disease ("candidate genes") to be associated with neurodevelopmental disorders (NDDs). Such variants are typically identified through genome wide screening approaches in individuals NDDs but without a clear diagnostic variant. AutoCaSc accounts for variant-specific parameters (conservation, "in silico" predictions), gene specific parameters (gene constraint, expression pattern, protein interactions), segregation of the variant and the overall interplay between these parameters.
            
            __What do the 4 subscores stand for?__  
            - __Variant Attributes (6 points max):__ These include conservation (GERP++), "in silico" predictions (MutationAssessor, MutationTaster, Sift), splice site predictions (MaxEntScan, AdaBoost, RandomForest) and expected impact (VEP).
            - __Gene Attributes (1 point max):__ These are gene constraint parameters from gnomAD; LOUEF for loss of function variants, Z for missense variants.
            - __Inheritance (2 points max):__ These points depend on inheritance of the variant of interest and segregation of the variant in the family.
            - __Gene Plausibility (6 points max):__ These points are calculated based on the gene's expression pattern, protein-protein interactions, animal model phenotypes, published articles on PubMed, de novo variants in the gene linked to NDD and other sources.

            __How can I enter multiple (compound heterogyous) variants?__  
            Just enter all your variants of interest by separating them by a comma. If "compound heterozygous" is selected, webAutoCaSc will automatically match variants in the same gene and process them as corresponding compound heterozygous variants.
            
            __What do the inheritance options stand for?__  
            - __De novo:__ De novo variants are identified only in the index and have not been inherited from the parents.
            - __Inherited dominant:__ In case of an inherited dominant variant, the variant of interested has been inherited by an equally affected parent.
            - __Homozygous recessive:__ Variant is identified in homozygous state in the index and in heterozygous state in both healthy parents.
            - __X-linked:__ X-linked variants are being inherited from the heterozygous carrier mother and cause a phenotype in a male descendant as he has only one affected allele and no healthy second allele to compensate. De novo variants on the X-chromosome, for both female and male index individuals, are account for in the de novo inheritance option.
            - __Compound heterozygous:__ Compound heterozygous variants are two different variants in the same gene but on different alleles. Each is inherited from only one heterozygous carrier parent.
            - __Unknown:__ The "unknown" option can be used if information on the parents and thus on segregation is missing.
            
            __What does _webAutoCaSc_ stand for?__  
            _AutoCaSc_ stands for __Auto__mated __Ca__ndidate __Sc__ore. We use the __web__ prefix to distinguish the command line interface (CLI) from the webapplication running the AutoCaSc script. The __Ca__ndidate __Sc__ore principle has been previously descrbed (Büttner et al. bioRxiv. 2019). 
            
            __Can webAutoCaSc be used for other phenotypes as well?__  
            AutoCaSc has been developed to work for NDDs. We don't recommend using it for other phenotypes. We are planning a generalized phenotype agnostic framework for future updates.
            """),
        ],
        # style={"min-height": "calc(99vh - 150px)"}
    ),
    # html.Div(style={"height": "1vh"}),
    footer
    ])

impressum_page = html.Div(
    [
        navbar,
        dbc.Container([
            dbc.Row([
                html.H2("Impressum"),
                dbc.Button("EN", id="impressum_language_button")
            ],
            justify="between",
            style={
                "padding-left": "18px",
                "padding-right": "18px",
            }),
            html.Br(),
            html.Div([
                dcc.Markdown("""
                            Gemäß § 28 BDSG widersprechen wir jeder kommerziellen Verwendung und Weitergabe der Daten.\n
                            __Verantwortunsbereich__:\n
                            Das Impressum gilt nur für die Internetpräsenz unter der Adresse: https://autocasc.uni-leipzig.de\n
                            __Abgrenzung__:\n
                            Die Web-Präsenz ist Teil des WWW und dementsprechend mit fremden, sich jederzeit wandeln könnenden Web-Sites verknüpft, die folglich auch nicht diesem Verantwortungsbereich unterliegen und für die nachfolgende Informationen nicht gelten. Dass die Links weder gegen Sitten noch Gesetze verstoßen, wurde genau ein Mal geprüft (bevor sie hier aufgenommen wurden).\n
                            __Diensteanbieter__:\n
                            Johann Lieberwirth und Rami Abou Jamra\n
                            __Ansprechpartner für die Webseite__:\n
                            Johann Lieberwirth (johann.lieberwirth@medizin.uni-leipzig.de)\n
                            __Verantwortlicher__:\n
                            Rami Abou Jamra (rami.aboujamra@medizin.uni-leipzig.de)\n
                            __Anschrift__:\n
                            Sekretariat\n
                            Philipp-Rosenthal-Str. 55\n
                            04103 Leipzig\n
                            Telefon: 0341 - 97 23800\n
                            __Urheberschutz und Nutzung__:\n
                            Der Urheber räumt Ihnen ganz konkret das Nutzungsrecht ein, sich eine private Kopie für persönliche Zwecke anzufertigen. Nicht berechtigt sind Sie dagegen, die Materialien zu verändern und /oder weiter zu geben oder gar selbst zu veröffentlichen.
                            Wenn nicht ausdrücklich anders vermerkt, liegen die Urheberrechte bei Johann Lieberwirth
                            Datenschutz Personenbezogene Daten werden nur mit Ihrem Wissen und Ihrer Einwilligung erhoben. Auf Antrag erhalten Sie unentgeltlich Auskunft zu den über Sie gespeicherten personenbezogenen Daten. Wenden Sie sich dazu bitte an den Administrator.\n
                            __Keine Haftung__:\n
                            Die Inhalte dieses Webprojektes wurden sorgfältig geprüft und nach bestem Wissen erstellt. Aber für die hier dargebotenen Informationen wird kein Anspruch auf Vollständigkeit, Aktualität, Qualität und Richtigkeit erhoben. Es kann keine Verantwortung für Schäden übernommen werden, die durch das Vertrauen auf die Inhalte dieser Website oder deren Gebrauch entstehen.\n
                            __Schutzrechtsverletzung__:\n
                            Falls Sie vermuten, dass von dieser Website aus eines Ihrer Schutzrechte verletzt wird, teilen Sie das bitte umgehend per elektronischer Post mit, damit zügig Abhilfe geschafft werden kann. Bitte nehmen Sie zur Kenntnis: Die zeitaufwändigere Einschaltung eines Anwaltes zur für den Diensteanbieter kostenpflichtigen Abmahnung entspricht nicht dessen wirklichen oder mutmaßlichen Willen.\n
                            \n
                            lt. Urteil vom 12. Mai 1998 - 312 O 85/98 - "Haftung für Links" hat das Landgericht Hamburg entschieden, dass man durch die Anbringung eines Links, die Inhalte der gelinkten Seite ggf. mit zu verantworten hat. Dies kann nur dadurch verhindert werden, dass man sich ausdrücklich von diesen Inhalten distanziert.
                            'Hiermit distanzieren wir uns ausdrücklich von allen Inhalten aller gelinkten Seiten auf unserer Website und machen uns diese Inhalte nicht zu eigen. Diese Erklärung gilt für alle auf unsere Website angebrachten Links.'
                            \n
                            © Copyright 2021
                """)
            ],
            id="impressum_text")
        ],
        # style={"min-height": "calc(99vh - 150px)"}
        ),
        # html.Div(style={"height": "1vh"}),
        footer
    ]
)



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
    html.Div(id="loading_output"),
    impressum_page,
])


########## FRONTEND ##########
# basic string_formatting and initiating AutoCaSc intsances in order to check if input is ok
def parse_input(input):
    return [input.split(",")[i].strip() for i in range(len(input.split(",")))]

def get_display_variant(_variant, n_chars=25):
    if len(_variant) < (n_chars + 5):
        return _variant
    else:
        return _variant[:n_chars] + "..."


def get_results_page(results_memory):
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
                        dbc.Row(
                            [
                                dbc.Tabs(
                                    tab_list,
                                    id="card_tabs",
                                    card=True,
                                    active_tab=initial_tab,
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
                        html.P(id="card_content"),
                        style={"padding-bottom": "0"}
                    )
                ])
            ],
                # style={"min-height": "calc(99vh - 135px)"}
        ),
            # html.Div(style={"height": "1vh"}),
            footer,
        ])
        return results_page


def input_ok(instances):
    for instance in instances:
        if instance.variant_format == "incorrect":
            return False
    return True


def get_error(error_code):
    error_dict = {
        201: ("Error: The reference sequence does not match GRCh37!", "warning"),
        301: ("Error: Could not identify the corresponding compound heterozygous variant!", "danger"),
        400: ("Error: Could not process variant. Please try VCF annotation or HGVS annotation using HGNC gene symbol!",
              "danger"),
        496: ("Error: The alternative sequence matches the GRCh37 reference sequence!", "danger"),
        498: ("Error: You have entered an intergenic variant.", "danger")
    }
    return error_dict.get(error_code) or (f"Some error occured. Code {error_code}", "warning")


def get_badge(status_code, i=None):
    if status_code == 200:
        return None
    else:
        error_message, error_color = get_error(status_code)
        return_badge = html.Div(
            [
                dbc.Badge(status_code, color=error_color, className="mx-2", id=f"return_badge_{i}"),
                dbc.Tooltip(error_message, target=f"return_badge_{i}")
            ]
        )
    return return_badge


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
    ctx = dash.callback_context

    if "pathname" in ctx.triggered[0]['prop_id']:
        if pathname == "/about":
            return about_page
        if pathname == "/impressum":
            return impressum_page
        if "/faq" in pathname:
            return faq_page
        if "/search" in pathname:
            if results_memory is None:
                print("search page")
                return search_page
            else:
                print("results page")
                return get_results_page(results_memory)
        else:
            return landing_page
    else:
        print("results page")
        return get_results_page(results_memory)


@app.callback(
    Output("variant_input", "valid"),
    Output("variant_input", "invalid"),
    Output("variant_queue_input", "data"),
    [Input("variant_input", "n_blur")],
    [State("variant_input", "value")]
)
def check_user_input(_, user_input):
    if not user_input == None:
        variants = parse_input(user_input)
        variant_instances = [AutoCaSc(_variant, mode="web") for _variant in variants]
        if input_ok(variant_instances):
            variant_queue = {}
            variant_queue["instances"] = {_instance.__dict__.get("variant"): _instance.__dict__ for _instance in variant_instances}
            return True, False, variant_queue
        else:
            return False, True, None
    raise PreventUpdate


@app.callback(
    Output("url", "pathname"),
    Input("search_button", "n_clicks"),
    State("variant_queue_input", "data"),
    State("variant_queue_url", "data"),
    State("inheritance_input", "value")
)
def search_button_click(n_clicks, variant_queue_input, variant_queue_url, inheritance):
    variant_queue = variant_queue_input or variant_queue_url
    if not n_clicks is None and not variant_queue is None and not inheritance is None:
        variants = [variant_queue.get("instances").get(_key).get("variant") for _key in variant_queue.get("instances").keys()]
        url_suffix = quote(f"/search/inheritance={inheritance}/variants={'%'.join(variants)}")
        return url_suffix
    else:
        raise PreventUpdate


@app.callback(
    Output("query_memory", "data"),
    Input("url", "pathname"),
    State("query_memory", "data")
)
def interpret_url_inheritance(pathname, query_memory):
    if "search" in pathname and query_memory is None:
        inheritance = unquote(pathname).split("inheritance=")[-1].split("/")[0]
        query_data = {}
        query_data["inheritance"] = inheritance
        return query_data
    raise PreventUpdate


@app.callback(
    Output("variant_queue_url", "data"),
    Input("url", "pathname"),
    State("variant_memory", "data")
)
def interpret_url_variants(pathname, variant_memory):
    if "search" in pathname and variant_memory is None:
        variants = unquote(pathname).split("variants=")[-1].split("%")
        variant_instances = [AutoCaSc(_variant, mode="web") for _variant in variants]
        variant_queue = {}
        variant_queue["instances"] = {_instance.__dict__.get("variant"): _instance.__dict__ for _instance in
                                      variant_instances}
        return variant_queue
    elif "search" in pathname:
        variants = unquote(pathname).split("variants=")[-1].split("%")
        if any([_variant not in variant_memory.get("instances").keys() for _variant in variants]):
            variant_instances = [AutoCaSc(_variant, mode="web") for _variant in variants]
            variant_queue = {}
            variant_queue["instances"] = {_instance.__dict__.get("variant"): _instance.__dict__ for _instance in
                                          variant_instances}
            return variant_queue
    raise PreventUpdate


def show_other_variant_column(results_memory):
    if any([results_memory.get("instances").get(_variant).get("inheritance") == "comphet"
            for _variant in results_memory.get("instances").keys()]):
        return True
    else:
        return False


@app.callback(
    Output("card_content", "children"),
    Input("card_tabs", "active_tab"),
    State("results_memory", "data"),
)
def get_tab_card(active_tab, results_memory):
    cell_style = {
        "padding": "5px",
        "padding-left": "12px"
    }
    if active_tab == "overview_tab":
        if show_other_variant_column(results_memory):
            other_variant_column_header = html.Th("Corresponding Variant")
        else:
            other_variant_column_header = None

        overview_table_header = html.Thead(
            html.Tr(
                [
                    html.Th("Variant"),
                    html.Th("Candidate Score"),
                    html.Th("HGVS"),
                    other_variant_column_header
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
                    else:
                        _hgvs += " "
                    _hgvs += _instance_attributes.get("hgvsp_change")
            else:
                _hgvs = None

            if show_other_variant_column(results_memory):
                other_variant_column = html.Th(
                    html.P(_instance_attributes.get("other_variant"),
                           id=f"other_variant_target_{i}"),
                           style=cell_style)
                tooltips.append(dbc.Tooltip(f"The corresponding compound heterozygous variant.",
                                            target=f"other_variant_target_{i}"))
            else:
                other_variant_column = None

            overview_table_rows.append(
                html.Tr(
                    [
                        html.Th(dbc.Row([dbc.Col(_variant,
                                                 width="auto"),
                                         dbc.Col(get_badge(_instance_attributes.get("status_code"), i),
                                                 width="auto")],
                                        justify="start",
                                        no_gutters=True),
                                style=cell_style),
                        html.Th(_instance_attributes.get("candidate_score_v1"),
                                style=cell_style),
                        html.Th(html.P(_hgvs,
                                       id=f"hgvs_target_{i}"),
                                style=cell_style),
                        other_variant_column
                    ]
                )
            )
            tooltips.append(dbc.Tooltip(f"Transcript: {_instance_attributes.get('transcript')}",
                            target=f"hgvs_target_{i}"))

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
                                style={"margin-bottom": "0",
                                       "padding-bottom": "0"}
                            ),
                            className="col-12 col-md-6"
                        )
                    ],
                style={"margin-bottom": "0",
                       "padding-bottom": "0"}
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
                dbc.Tooltip(f"{_instance_attributes.get('explanation_dict').get('impact')}, "
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
                dbc.Tooltip("Points attributed for data in literature databases, animal models, expression pattern, "
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
                # html.Br(),
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
        else:
            error_message, error_color = get_error(status_code)
            return dbc.Alert(error_message, color=error_color)


@app.callback(
    Output("about_text", "children"),
    Output("about_language_button", "children"),
    Input("about_language_button", "n_clicks"),
    State("about_language_button", "children")
)
def get_about_text(n_clicks, language):
    if not n_clicks:
        raise PreventUpdate
    if language == "DE":
        return [
            html.P("AutoCaSc ist ein Skript zum automatisierten Bewerten von Kandidatenvarianten in Fällen neuronaler Entwicklungsverzögerung. Es ist ausschließlich für Forschungszwecke zu benutzen. Die Annotation der Varianten erfolgt mit der REST API von VEP (ensembl, (1)). Zur Berechnung der Kandidatenpunktzahl (Candidate score, CaSc) werden 12 verschiedene Parameter einbezogen. Dies sind die Art der Vererbung (z.B. de novo), Genattribute wie pLI & Z-Score (gnomAD (2)), Expressionsmuster (GTEx (3)), in silico Analysen, Proteininteratkionsdatnebanken (StringDB (4)), Tierdatenbanken (MGI (5)), Literaturdatenbanken (Pubtator Central (6)),  sowie weitere (PsymuKB (7), DisGeNET (8)). Die maximal erreichbare Punktzahl sind 15 Punkte. Je höher der erreichte Punktwert, desto plausibler scheint die aus den zugrundeliegenden Daten errechnete Pathogenität der Variante mit Blick auf neuronale Entwicklungsverzögerung.",
                   style={"text-align": "justify"}),
            html.Br(),
            html.P("Bei Fragen und Anmerkungen kontaktieren Sie bitte johann.lieberwirth@medizin.uni-leipzig.de oder rami.aboujamra@medizin.uni-leipzig.de.",
                   style={"text-align": "justify"})
               ], "EN"
    else:
        return [
                html.P("AutoCaSc is a tool for quantifying the plausibility of candidate variants for Neurodevelopmental Disorders (NDD). AutoCaSc is intended to be used on impactful rare variants in a research setting. In its current version, 12 parameters are counted in, achieving a maximum of 15 points. User inputs are the identified variant in a standard HGVS/VCF format together with segregation aspects (de novo, recessive, dominant and X-chromosomal). We use the Ensembl VEP REST API (1) to annotate variant attributes (e.g. variant consequence, allele frequency from gnomAD, in silico predictions) and gene based scores dependent on inheritance mode (e.g. high Z-score is of relevant for de novo missense) from dbNSFP. Other attributes were previously labor intensive and predisposed to variability. These included important categories like expression in the nervous system, neuronal functions, co-expression and protein interactions, search for relevant literature, model organisms and observations in screening studies. As an objective approach we now searched a defined set of databases (gnomAD (2), GTEx (3), STRING (4), MGI (5), PubTator (6), PsyMuKB (7), DisGeNET (8)) and generated empirical cut-offs for each category by comparing the respective readout between a manually curated list of known NDD genes from the SysID database (9) and a list of genes not involved in NDD.",
                       style={"text-align": "justify"}),
                html.Br(),
                html.P("Feel free to contact johann.lieberwirth@medizin.uni-leipzig.de or rami.aboujamra@medizin.uni-leipzig.de in case you have further questions or in case you have found a bug.",
                       style={"text-align": "justify"})
            ], "DE"


@app.callback(
    Output("impressum_text", "children"),
    Output("impressum_language_button", "children"),
    Input("impressum_language_button", "n_clicks"),
    State("impressum_language_button", "children")
)
def get_impressum_text(n_clicks, language):
    if not n_clicks:
        raise PreventUpdate
    if language == "DE":
        return dcc.Markdown("""
                            Gemäß § 28 BDSG widersprechen wir jeder kommerziellen Verwendung und Weitergabe der Daten.\n
                            __Verantwortunsbereich__:\n
                            Das Impressum gilt nur für die Internetpräsenz unter der Adresse: https://autocasc.uni-leipzig.de\n
                            __Abgrenzung__:\n
                            Die Web-Präsenz ist Teil des WWW und dementsprechend mit fremden, sich jederzeit wandeln könnenden Web-Sites verknüpft, die folglich auch nicht diesem Verantwortungsbereich unterliegen und für die nachfolgende Informationen nicht gelten. Dass die Links weder gegen Sitten noch Gesetze verstoßen, wurde genau ein Mal geprüft (bevor sie hier aufgenommen wurden).\n
                            __Diensteanbieter__:\n
                            Johann Lieberwirth und Rami Abou Jamra\n
                            __Ansprechpartner für die Webseite__:\n
                            Johann Lieberwirth (johann.lieberwirth@medizin.uni-leipzig.de)\n
                            __Verantwortlicher__:\n
                            Rami Abou Jamra (rami.aboujamra@medizin.uni-leipzig.de)\n
                            __Anschrift__:\n
                            Sekretariat\n
                            Philipp-Rosenthal-Str. 55\n
                            04103 Leipzig\n
                            Telefon: 0341 - 97 23800\n
                            __Urheberschutz und Nutzung__:\n
                            Der Urheber räumt Ihnen ganz konkret das Nutzungsrecht ein, sich eine private Kopie für persönliche Zwecke anzufertigen. Nicht berechtigt sind Sie dagegen, die Materialien zu verändern und /oder weiter zu geben oder gar selbst zu veröffentlichen.
                            Wenn nicht ausdrücklich anders vermerkt, liegen die Urheberrechte bei Johann Lieberwirth
                            Datenschutz Personenbezogene Daten werden nur mit Ihrem Wissen und Ihrer Einwilligung erhoben. Auf Antrag erhalten Sie unentgeltlich Auskunft zu den über Sie gespeicherten personenbezogenen Daten. Wenden Sie sich dazu bitte an den Administrator.\n
                            __Keine Haftung__:\n
                            Die Inhalte dieses Webprojektes wurden sorgfältig geprüft und nach bestem Wissen erstellt. Aber für die hier dargebotenen Informationen wird kein Anspruch auf Vollständigkeit, Aktualität, Qualität und Richtigkeit erhoben. Es kann keine Verantwortung für Schäden übernommen werden, die durch das Vertrauen auf die Inhalte dieser Website oder deren Gebrauch entstehen.\n
                            __Schutzrechtsverletzung__:\n
                            Falls Sie vermuten, dass von dieser Website aus eines Ihrer Schutzrechte verletzt wird, teilen Sie das bitte umgehend per elektronischer Post mit, damit zügig Abhilfe geschafft werden kann. Bitte nehmen Sie zur Kenntnis: Die zeitaufwändigere Einschaltung eines Anwaltes zur für den Diensteanbieter kostenpflichtigen Abmahnung entspricht nicht dessen wirklichen oder mutmaßlichen Willen.\n
                            \n
                            lt. Urteil vom 12. Mai 1998 - 312 O 85/98 - "Haftung für Links" hat das Landgericht Hamburg entschieden, dass man durch die Anbringung eines Links, die Inhalte der gelinkten Seite ggf. mit zu verantworten hat. Dies kann nur dadurch verhindert werden, dass man sich ausdrücklich von diesen Inhalten distanziert.
                            'Hiermit distanzieren wir uns ausdrücklich von allen Inhalten aller gelinkten Seiten auf unserer Website und machen uns diese Inhalte nicht zu eigen. Diese Erklärung gilt für alle auf unsere Website angebrachten Links.'
                            \n
                            © Copyright 2021
                """), "EN"
    else:
        return [
                html.P("The Institute for Human Genetics (University Medical Center Leipzig) makes no representation about the suitability or accuracy of this software or data for any purpose, and makes no warranties, including fitness for a particular purpose or that the use of this software will not infringe any third party patents, copyrights, trademarks or other rights."),
                html.Br(),
                dcc.Markdown("""
                __Responsible for this website__:\n
                Johann Lieberwirth (johann.lieberwirth@medizin.uni-leipzig.de)\n
                __Responsible for this project__:\n
                Rami Abou Jamra (rami.aboujamra@medizin.uni-leipzig.de)\n
                __Address__:\n
                Sekretariat\n
                Philipp-Rosenthal-Str. 55\n
                04103 Leipzig\n
                GERMANY\n
                Telefon: 0341 - 97 23800\n""")
            ], "DE"

########## BACKEND ##########
def score_variants(instances, inheritance):
    instances_processed = []
    if inheritance == "comphet":
        variant_gene_df = pd.DataFrame()
        for i, _instance in enumerate(instances):
            _instance.inheritance = inheritance
            if _instance.__dict__.get("status_code") == 200:
                variant_gene_df.loc[i, "variant"] = _instance.__dict__.get("variant")
                variant_gene_df.loc[i, "gene_symbol"] = _instance.__dict__.get("gene_symbol")
                variant_gene_df.loc[i, "instance"] = _instance
            else:
                instances_processed.append(_instance)
        for _gene in variant_gene_df.gene_symbol.unique():
            df_chunk = variant_gene_df.loc[variant_gene_df.gene_symbol == _gene].reset_index(drop=True)
            if len(df_chunk) == 2:
                for i in range(2):
                    _instance = copy.deepcopy(df_chunk.loc[abs(1-i), "instance"])
                    _instance.other_autocasc_obj = df_chunk.loc[abs(0-i), "instance"]
                    _instance.calculate_candidate_score()
                    instances_processed.append(_instance)
            else:
                for i in range(len(df_chunk)):
                    lonely_instance = df_chunk.loc[i, "instance"]
                    lonely_instance.__dict__["status_code"] = 301
                    instances_processed.append(lonely_instance)
                    # raise IOError(f"combination of comphet variants unclear for {_gene}")
    else:
        for _instance in instances:
            _instance.inheritance = inheritance
            if _instance.__dict__.get("status_code") == 200:
                _instance.calculate_candidate_score()
            instances_processed.append(_instance)
    return instances_processed


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
    Input("variant_queue_input", "data"),
    Input("variant_queue_url", "data"),
    [State("variant_memory", "data")]
)
def retrieve_variant_data(variant_queue_input, variant_queue_url, variant_memory):
    variant_queue = variant_queue_input or variant_queue_url
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
    Input("query_memory", "data"),
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
    # todo check comphet behavior
    if not n_cklicks:
        raise PreventUpdate
    df = pd.DataFrame()
    for i, _variant in enumerate(results_memory.get("instances").keys()):
        _instance = results_memory.get("instances").get(_variant)
        try:
            if _instance.get("vcf_string") in df["vcf_format_2"].to_list():
                continue
        except KeyError:
            pass
        if _instance.get("inheritance") == "comphet":
            comphet = True
            _other_instance = _instance.get("other_autocasc_obj")
        else:
            comphet = False
            _other_instance = {}
        df.loc[i, "hgnc_symbol"] = _instance.get("gene_symbol")
        df.loc[i, "transcript"] = _instance.get("transcript")
        df.loc[i, "vcf_format_1"] = _instance.get("vcf_string")
        df.loc[i, "vcf_format_2"] = _other_instance.get("vcf_string")
        df.loc[i, "cDNA_1"] = _instance.get("hgvsc_change")
        df.loc[i, "cDNA_2"] = _other_instance.get("hgvsc_change")
        df.loc[i, "amino_acid_1"] = _instance.get("hgvsp_change")
        df.loc[i, "amino_acid_2"] = _other_instance.get("hgvsp_change")
        df.loc[i, "var_1_full_name"] = f"{_instance.get('transcript')}:{_instance.get('hgvsc_change')} {_instance.get('hgvsp_change')}"
        if comphet:
            df.loc[i, "var_2_full_name"] = f"{_other_instance.get('transcript')}:{_other_instance.get('hgvsc_change')} {_other_instance.get('hgvsp_change')}"
        else:
            df.loc[i, "var_2_full_name"] = ""
        df.loc[i, "inheritance"] = _instance.get("inheritance")
        df.loc[i, "candidate_score"] = _instance.get("candidate_score_v1")
        df.loc[i, "literature_plausibility"] = _instance.get("literature_score")
        df.loc[i, "inheritance_score"] = _instance.get("inheritance_score")
        if comphet:
            df.loc[i, "variant_attribute_score"] = round(mean([_instance.get("variant_score"),
                                                               _other_instance.get("variant_score")]), 2)
        else:
            df.loc[i, "variant_attribute_score"] = _instance.get("variant_score")
        df.loc[i, "gene_attribute_score"] = _instance.get("gene_attribute_score")

    data = io.StringIO()
    df.to_csv(data, sep="\t", decimal=",")
    data.seek(0)
    return dict(content=data.getvalue(), filename="AutoCaSc_results.tsv")


if __name__ == '__main__':
    app.run_server(debug=True,
                   dev_tools_hot_reload=False,
                   host='0.0.0.0',
                   port=5000)