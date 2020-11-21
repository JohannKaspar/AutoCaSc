import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from statistics import mean
from dash.dependencies import Input, Output, State
from GTEx import GeneExpression  # used for expression visualization
from gnomAD import GnomADQuery
import dash_table
import re
# from dash.exceptions import PreventUpdate
from AutoCaSc import AutoCaSc, AUTOCASC_VERSION
from urllib.parse import unquote, quote
import pandas as pd
from tools import filterTheDict

DASHUI_VERSION = "w0.02"

# CSS stylsheets sind im Ordner "assets"
app = dash.Dash(__name__)

###### (STYLE)DICTIONARIES ######
colors = {
    'background': '#ffffff',
    'accent': '#7CD3CB',
    "text": "#111111",
    "brain": "#fadf63",
    "other": "#3877B6",
    "dark_grey": "#5F7470",
    "light_grey": "#889696",
    "red": "#E15554",
    "green": "#79B473",
    "yellow": "#F2C57C",
    "doughnut_blue": "#5f80b9",
    "light_blue": "#cad5e8",
    "doughnut_red": "#a8504e",
    "light_red": "#e7cbcb",
    "doughnut_orange": "#e0934d",
    "light_orange": "#f4d8be",
    "doughnut_yellow": "#efbf2c",
    "light_yellow": "#fae9b8",
}

# list of parameters to show in variant results table from variant_annotation_result
variant_parameters_to_show = [
    "id",
    "gene_symbol",
    "vcf_string",
    "transcript",
    "hgvsc_change",
    "gene_id",
    "hgvsp_change",
    "consequence",
    "impact",
    "oe_lof_interval",
    "mis_z",
    "oe_mis_interval",
    "pLI",
    "gnomad_frequency",
    "MAF"
]
# ToDo beide dictionarys mergen, maxentscan_consequence/ada_consequence/rf_consequence mit anzeigen?

# dict of labels for variant results table
variant_parameters_labels = {
    "gene_symbol": "Gen",
    "id": "Variant",
    "vcf_string": "VCF",
    "transcript": "Transcript",
    "hgvsp_change": "HGVSP",
    "hgvsc_change": "HGVSC",
    "gene_id": "gencodeID",
    "impact": "Impact",
    "consequence": "Consequence",
    "oe_lof_interval": "oe LoF",
    "mis_z": "Z",
    "oe_mis_interval": "oe mis",
    "pLI": "pLI",
    "gnomad_frequency": "gnomAD Ex Freq",
    "MAF": "gnomAD MAF",
    "hgnc_id": "HGNC ID",
    "omim_id": "OMIM ID"
}

# list of parameters to show in prediction results table from variant_annotation_result
prediction_scores_to_show = [
    "sift_converted_rankscore",
    "mutationtaster_converted_rankscore",
    "mutationassessor_rankscore",
    "cadd_phred",
    "gerp_rs_rankscore",
    "polyphen_prediction",
    "ada_score",
    "rf_score",
    "maxentscan_decrease"
]

# dict of labels for prediction results table
prediction_scores_labels = {
    "cadd_phred": "CADD (phred)",
    "mutationtaster_converted_rankscore": "MutationTaster RS",
    "polyphen_prediction": "Polyphen prediction",
    "gerp_rs_rankscore": "GERP++ RS",
    "mutationassessor_rankscore": "MutationAssessor RS",
    "sift_converted_rankscore": "Sift RS",
    "ada_score": "Splicing: Ada Boost",
    "rf_score": "Splicing: Random Forest",
    "maxentscan_decrease": "Splicing: MaxEntScan Decrease"
}


###### HELPER FUNCTIONS ######
# function that looks for element in a given list by a given index, if nonexistent returns "None"
def safe_get(l, idx, default=None):
    try:
        return l[idx]
    except IndexError:
        return default
    except AttributeError:
        return default


# creates static half of left column
def generate_description_card():
    return html.Div(
        id="description-card",
        children=[
            html.H5("AutoCaSc"),
            html.H3("Welcome to the Automated Candidate Scoring Tool"),
            html.Div(
                id="intro",
                children="Score your candidate genes by entering a Variant in VCF format. The result is calculated based on aggregation of variant prediction scores and gene-specific parameters.",
            ),
            html.Hr()
        ],
    )


# creates dynamic half of left column with user inputs
def generate_control_card(inheritance_dd=None,
                          family_history_radio=[]):
    return html.Div(
        id="control-card",
        children=[
            html.P("Enter a variant"),
            dcc.Textarea(
                id='text_input',
                draggable=False,
                placeholder='Enter a variant in VCF format...',
                style={"width": "100%",
                       "max-width": "100%",
                       "min-width": "100%"},
            ),
            html.Hr(),

            # outputs examples for possible inputs
            html.P("Examples"),
            html.Div(
                id="example_card",
                children=[
                    dcc.Link("ChrX-153040798-C-A", href="/X:153040798:C:A", className="example_link", id="example_a"),
                    dcc.Link("1:7725246:G:A", href="/1:7725246:G:A", className="example_link", id="example_b"),
                    dcc.Link("10 123809984 C CCCTC", href="/10:123809984:C:CCCTC", className="example_link",
                             id="example_c"),
                    dcc.Link("19 53959842 xxx -", href="/19:53959842:xxx:-", className="example_link", id="example_d"),
                    dcc.Link("AGT:c.803T>C", href="/AGT:c.803T>C", className="example_link", id="example_e"),
                    dcc.Link("ENST00000003084: c.1431_1433delTTC", href="/ENST00000003084:c.1431_1433delTTC",
                             className="example_link", id="example_f"),
                ]
            ),
            html.Hr(),
            # html.Br(),

            # dropdown menu for reference genome
            html.P("Select reference genome"),
            dcc.Dropdown(
                id="dd_reference_genome",
                options=[{"label": "GRCh37", "value": "GRCh37"}, {"label": "GRCh38", "value": "GRCh38"}],
                value="GRCh37",
            ),
            html.Hr(),

            # dropdown for initial inheritance options
            html.P("Select Inheritance"),
            dcc.Dropdown(
                id="inheritance_dropdown",
                options=[{"label": "de novo", "value": "de_novo"},
                         {"label": "homo (auto/recessive)", "value": "homo"},
                         {"label": "compound heterozygous", "value": "comphet"},
                         {"label": "X-linked", "value": "x_linked"},
                         {"label": "inherited autosomal dominant", "value": "ad_inherited"},
                         {"label": "other", "value": "other"}],
                value="other"
            ),
            html.Hr(style={"margin-bottom": "0.5rem"}),

            # choices for further specification of inheritance, depending on input from innheritance dropdown menu
            html.Div(
                id="family_history_radio_card",
                children=[
                    dcc.RadioItems(
                        id="family_history_radio"
                    ),
                    dcc.RadioItems(
                        id="other_impact"
                    )
                ]
            ),

            # container for search button and loading sign
            html.Div(
                id="search_loading_card",
                children=[
                    html.Div(
                        id="search_button_card",
                        children=html.Button(
                            id="search_button",
                            children="Start Search",
                            n_clicks=0,
                        ),
                        style={"width": "66%"}
                    ),
                    html.Div(
                        id="loading_spacer",
                        children=[html.Div(
                            id="loading_sign_card",
                            children=dcc.Loading(
                                id="loading_sign",
                                type="circle",
                                loading_state={"is_loading": True},
                                children=html.Div(id="loading-output-1"),
                                className="container"
                            ),
                        )
                        ],
                        className="four columns",
                    )
                ]),

            # version number + copyrights
            html.Small(f"version {AUTOCASC_VERSION}.{DASHUI_VERSION}", id="version_container")
            # html.Div(
            #     id="version_container",
            #     children=[
            #
            #     ]
            # )
        ],
        # style={"width":"100%"},
    )


# clears page and outputs a message in case of an error
def clear_page(status_code, variant_input=""):
    error_dict = {
        400: ("Error 400: Variant incorrect. Please try VCF annotation or HGVS annotation using HGNC gene symbol!",
              colors["red"]),
        401: ("Error 401: Variant incorrect. Please try VCF annotation or HGVS annotation using HGNC gene symbol!",
              colors["red"]),  # regex doesn't match
        402: ("Error 402: Can't score a compound heterozygous variant without the correspoing variant!",
              colors["red"]),
        497: ("Error 497: Variant incorrect. Please try VCF annotation or HGVS annotation using HGNC gene symbol!",
              colors["red"]),
        498: ("You have entered an intergenic variant.", colors["red"]),
        496: ("The alternative sequence matches the reference sequence!", colors["red"]),
        300: (" ", colors["light_grey"])
    }

    error_message, error_color = error_dict.get(status_code) or (
    f"There has been a {status_code} error.", colors["red"])

    return ["",
            variant_input,
            {},
            {},
            [
                # status bar
                html.Div(
                    id="status_bar",
                    children=error_message,
                    style={"backgroundColor": error_color},
                ),

                # CaSc calculation results card
                html.Div(
                    id="casc_results_card",
                    style={"borderColor": "#f1f1f1"},
                    hidden=True,
                    children=[
                        html.Div(
                            html.B("Candidate Score", className="card_title"),
                            className="card_header"
                        ),
                        html.Div(
                            id="casc_results",
                            children=[
                                get_table("calculation_results_table", {"": ""}, {"": ""}, {"": ""})
                            ]
                        ),
                        html.Div(
                            id="copy_results_field_container"
                        )
                    ],
                ),

                # variant annotation results
                html.Div(
                    html.Div(
                        id="variant_card",
                        className="pretty_container",
                        hidden=True,
                        children=[
                            html.Div([
                                html.B("Variant Annotation", className="card_title"),
                            ],
                                className="card_header"
                            ),
                            html.Div(id="variant_annotation",
                                     children=[
                                         get_table("empty_variant", {"": ""}, {"": ""}, {"": ""})
                                     ]
                                     ),
                            # html.Hr(),
                        ],
                    ),
                ),

                # tissue expression graph
                html.Div(
                    id="tissue_expression_card",
                    className="pretty_container",
                    hidden=True,
                    children=[
                        html.Div([
                            html.B("Tissue Expression", className="card_title"),
                        ],
                            className="card_header"
                        ),
                        html.Hr(),
                        dcc.Graph(id="tissue_expression_graph", config={'displayModeBar': False}),
                    ],
                ),
            ]
            ]


def get_status_bar(status_code):
    if status_code == 200:
        status_bar_color = {
            "backgroundColor": colors["green"]
        }
        status_bar_text = "Annotation successful!"
        casc_card_border_color = {
            "borderColor": colors["green"]
        }
    elif status_code == 201:
        status_bar_color = {
            "backgroundColor": colors["yellow"]
        }
        status_bar_text = "CAVE! Reference Sequences don't match!"
        casc_card_border_color = {
            "borderColor": colors["yellow"]
        }

    return status_bar_color, status_bar_text, casc_card_border_color


# function that takes in annotation results, filters relevant values and returns a formattet table (dcc.DataTable)
def get_table(table_name, data, param_list=[], label_dict={}, orientation="vertical", editable=False):
    # in case of no specification use all given values for table creation
    if data is not None and data != {'': ''}:
        if param_list == []:
            param_list = list(data.keys())
        if label_dict == {}:
            for param in param_list:
                label_dict[param] = param

    if orientation == "horizontal":
        table_columns = []
        table_data = {}

        for k in param_list:

            table_columns.append(
                {"name": label_dict.get(k), "id": k}
            )

            if data.get(k) is not None:
                table_data[k] = data.get(k)
            else:
                table_data[k] = "-"
        table_data = [table_data]

    elif orientation == "vertical":
        table_data = []
        table_columns = [{"name": "Parameter", "id": "parameter"}, {"name": "Value", "id": "value"}]
        for k in param_list:
            if data.get(k) is not None and "nan" not in str(data.get(k)):
                table_data.append({"parameter": label_dict.get(k), "value": data.get(k)})
            else:
                table_data.append({"parameter": label_dict.get(k), "value": "-"})

    table = dash_table.DataTable(
        id=table_name,
        columns=table_columns,
        data=table_data,

        editable=editable,
        style_as_list_view=True,

        style_data_conditional=
        [
            {
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'
            }
        ],

        style_data=
        {
            'fontFamily': 'Acumin, "Helvetica Neue", sans-serif',
            'fontWeight': 'regular',
            'fontSize': '1.3rem'
        },

        style_header=
        {
            'backgroundColor': 'rgb(230, 230, 230)',
            'fontWeight': 'bold',
            "fontFamily": 'Acumin, "Helvetica Neue", sans-serif'
        },

        style_cell={
            'textAlign': "left",
            "padding": "5px"
        },

        fill_width=True,
        style_table=
        {
            'overflowX': 'scroll',
        }
    )
    return table


def get_casc_card(casc_table, explanation_button):
    return [casc_table,
            html.Br(),
            html.Div(
                [
                    explanation_button,
                    html.Button(["copy results"],
                                id="copy_button",
                                className="six columns"),
                ]
            )]


def get_casc_table(table_name, data, param_list=[], label_dict={}, orientation="vertical", editable=True):
    # in case of no specification use all given values for table creation
    if data:
        if param_list == []:
            param_list = list(data.keys())
        if label_dict == {}:
            for param in param_list:
                label_dict[param] = param

        table_data = []
        table_columns = [{"name": "Variant", "id": "parameter"},
                         {"name": str(data.get("candidate_score")), "id": "value"}]

        for k in param_list:
            if data.get(k) is not None:
                table_data.append({"parameter": label_dict.get(k).replace("_", " ").title(), "value": data.get(k)})
            else:
                table_data.append({"parameter": label_dict.get(k).replace("_", " ").title(), "value": "-"})

    table = dash_table.DataTable(
        id=table_name,
        columns=table_columns,
        data=table_data,

        editable=editable,
        style_as_list_view=True,

        style_data_conditional=
        [
            {
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'
            },
            # make first row with final result bold
            # {
            #     "if": {"row_index": 0},
            #     'backgroundColor': 'rgb(230, 230, 230)',
            #     'fontWeight': 'bold'
            # },
            {
                "if": {"row_index": 0},
                "backgroundColor": colors["light_blue"]
            },
            {
                "if": {"row_index": 1},
                "backgroundColor": colors["light_red"]
            },
            {
                "if": {"row_index": 2},
                "backgroundColor": colors["light_yellow"]
            },
            {
                "if": {"row_index": 3},
                "backgroundColor": colors["light_orange"]
            },

        ],

        style_data=
        {
            'fontFamily': 'Acumin, "Helvetica Neue", sans-serif',
            'fontWeight': 'regular',
            'fontSize': '1.5rem'
        },

        style_header=
        {
            'backgroundColor': 'rgb(230, 230, 230)',
            'fontWeight': 'bold',
            "fontFamily": 'Acumin, "Helvetica Neue", sans-serif',
            "visibile": "False"
        },

        style_cell={
            'padding': '5px',
            'textAlign': 'left'
        },

        fill_width=True,

    )
    return table


def get_casc_batch_table(variant_results):
    casc_results_dict = {}

    table_data = []
    table_columns = [{"name": "Variant", "id": "variant", "type": "text", "presentation": "markdown"},
                     {"name": "Candidate Score", "id": "candidate_score"},
                     {"name": "Inheritance Score", "id": "inheritance_score"},
                     {"name": "Gene Attribute Score", "id": "gene_attribute_score"},
                     {"name": "Variant Score", "id": "variant_score"},
                     {"name": "Literature Score", "id": "literature_score"}]

    for variant in variant_results.keys():
        casc_results_dict[variant] = variant_results.get(variant).candidate_score
    sorted_variants = [tup[0] for tup in sorted(casc_results_dict.items(), key=lambda item: item[1], reverse=True)]

    for variant in sorted_variants:
        if variant_results.get(variant) is not None:
            # link_label = re.sub(r"\s+", "", variant)
            table_data.append({"variant": f"[{variant}](http://autocasc.pythonanywhere.com/{quote(variant)})",
                               "candidate_score": variant_results.get(variant).candidate_score,
                               "inheritance_score": variant_results.get(variant).inheritance_score,
                               "gene_attribute_score": variant_results.get(variant).gene_attribute_score,
                               "variant_score": variant_results.get(variant).variant_score,
                               "literature_score": variant_results.get(variant).literature_score})
        else:
            table_data.append({"variant": "-",
                               "candidate_score": "-",
                               "inheritance_score": "-",
                               "gene_attribute_score": "-",
                               "variant_score": "-",
                               "literature_score": "-"})
        table = dash_table.DataTable(
            id="calculation_results_table",
            editable=False,
            columns=table_columns,
            data=table_data,
            style_as_list_view=True,
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                },

                {
                    "if": {"column_id": "inheritance_score"},
                    "backgroundColor": colors["light_blue"]
                },
                {
                    "if": {"column_id": "gene_attribute_score"},
                    "backgroundColor": colors["light_red"]
                },
                {
                    "if": {"column_id": "variant_score"},
                    "backgroundColor": colors["light_yellow"]
                },
                {
                    "if": {"column_id": "literature_score"},
                    "backgroundColor": colors["light_orange"]
                },

            ],

            style_data=
            {
                'fontFamily': 'Acumin, "Helvetica Neue", sans-serif',
                'fontWeight': 'regular',
                'fontSize': '1.5rem',
            },

            style_header=
            {
                'backgroundColor': 'rgb(230, 230, 230)',
                'fontWeight': 'bold',
                "fontFamily": 'Acumin, "Helvetica Neue", sans-serif'
            },

            style_cell={
                'padding': '5px',
                'textAlign': 'left'
            },
            fill_width=True,
        )
    return table


def get_explanation_button(instance):
    explanation_dict = instance.explanation_dict
    button = dcc.ConfirmDialogProvider(
        children=html.Button("Explanation",
                             className="six columns"),
        id="explanation_button",
        message=f"""
        Inheritance
        {explanation_dict["inheritance"]}
        ---------------------------------------------
        Gene attributes
        pLI & Z:              {explanation_dict["pli_z"]}
        ---------------------------------------------
        Variant attributes
        impact:              {explanation_dict["impact"]}
        in silico:            {explanation_dict["in_silico"]}
        conservation:   {explanation_dict["conservation"]}
        frequency:        {explanation_dict["frequency"]}
        ---------------------------------------------
        Literature research
        weighted sum:       {instance.literature_score}
        pubtator:           {instance.pubtator_score}
        gtex:               {instance.gtex_score}
        denovo rank score:  {instance.denovo_rank_score}
        disgenet:           {instance.disgenet_score}
        MGI:                {instance.mgi_score}
        StringDB:           {instance.string_score}
        """
    )
    return button


def get_clipboard_format_single(annotation_results):
    value_list = []
    key_list = [
        "gene_symbol",
        "candidate_score",
        "inheritance_score",
        "gene_attribute_score",
        "variant_score",
        "literature_score",
        "mgi_score",
        "pubtator_score",
        "gtex_score",
        "denovo_rank_score",
        "disgenet_score",
        "string_score"
    ]

    for key in key_list:
        value = annotation_results.get(key)
        if value is not None:
            value_list.append(str(value))
        else:
            value_list.append("")

    value_list.append("v" + str(annotation_results.get("version")).replace(".", "_"))

    return_string = "\t".join(value_list)
    return return_string


# def get_clipboard_format_batch(annotation_results):
#
#     value_list = []
#     key_list = [
#         "gene_symbol",
#         "candidate_score",
#         "inheritance_score",
#         "gene_attribute_score",
#         "variant_score",
#         "literature_score",
#         "mgi_score",
#         "pubtator_score",
#         "gtex_score",
#         "denovo_rank_score",
#         "disgenet_score",
#         "string_score"
#     ]
#
#     for key in key_list:
#         value = annotation_results.get(key)
#         if value is not None:
#             value_list.append(str(value))
#         else:
#             value_list.append("")
#
#     value_list.append("v" + str(annotation_results.get("version")).replace(".", "_"))
#
#     return_string = "\t".join(value_list)
#     return return_string

def get_output_elements(instance):
    # updated_variant_table = get_table("variant_annotation_table",
    #                                   instance.__dict__,
    #                                   variant_parameters_to_show,
    #                                   variant_parameters_labels)
    #
    # updated_prediction_table = get_table("prediction_score_table",
    #                                      instance.__dict__,
    #                                      prediction_scores_to_show,
    #                                      prediction_scores_labels)
    variant_parameters_labels.update(prediction_scores_labels)
    updated_variant_table = get_table("variant_annotation_table",
                                      instance.__dict__,
                                      variant_parameters_to_show + prediction_scores_to_show,
                                      variant_parameters_labels)

    updated_prediction_table = get_table("empty", {"": ""}, {"": ""}, {"": ""})
    clipboard_data = {"data": get_clipboard_format_single(instance.__dict__)}

    casc_results_children = get_casc_card(get_casc_table("calculation_results_table",
                                                         instance.__dict__,
                                                         ["inheritance_score", "gene_attribute_score", "variant_score",
                                                          "literature_score"],
                                                         editable=True),
                                          get_explanation_button(instance))

    expression_by_tissue, gtex_status_code = GeneExpression(instance.gene_id).get_expression()

    if gtex_status_code != 200:
        # in case of errors during GTEx data retrieval returns an empty graph
        updated_graph = {"data": [], "layout": []}
        # ToDo hier auch sagen was das Problem ist

    else:
        # creates graph with retrieved expression data
        updated_graph = get_boxplot(expression_by_tissue)

    return updated_variant_table, updated_prediction_table, casc_results_children, clipboard_data, updated_graph, gtex_status_code


def get_boxplot(expression_by_tissue):
    brain_x = [x for x in expression_by_tissue.keys() if "Brain" in x]
    brain_y = [round(expression_by_tissue[x], 3) for x in brain_x]
    other_x = [x for x in expression_by_tissue.keys() if x not in brain_x]
    other_y = [round(expression_by_tissue[x], 3) for x in other_x]

    # formatting labels
    try:
        brain_x = [re.sub(r" - ", " ", x.split("Brain - ")[1]) for x in brain_x]
        other_x = [re.sub(r" - ", " ", x) for x in other_x]
    except IndexError:  # this is due to different formats depending on source (file or gtex-api)
        brain_x = [re.sub(r"_", " ", x.split("Brain_")[1]) for x in brain_x]
        other_x = [re.sub(r"_", " ", x) for x in other_x]

    # positioning the plots
    fig = go.Figure()

    ## all plots are being added
    fig.add_trace(go.Box(
        y=brain_y,
        text=brain_x,
        hovertemplate="<b>%{text}</b><br>%{y}",
        hoverlabel={
            "bgcolor": colors["brain"]},
        name="Brain",
        quartilemethod="linear",
        line_color=colors["brain"],
        boxpoints="all",
        pointpos=0,
        marker={
            "color": colors["light_grey"],
            "size": 6,
            "line": {
                "color": colors["dark_grey"],
                "width": 1
            }
        },
    ))

    fig.add_trace(go.Box(
        y=other_y,
        text=other_x,
        hovertemplate="<b>%{text}</b><br>%{y}",
        hoverlabel={
            "bgcolor": colors["other"]},
        name="Other tissues",
        quartilemethod="linear",
        line_color=colors["other"],
        boxpoints="all",
        pointpos=0,
        marker={
            "color": colors["light_grey"],
            "size": 6,
            "line": {
                "color": colors["dark_grey"],
                "width": 1
            }
        },
    ))

    fig.update_layout(
        plot_bgcolor=colors["background"],
        height=450,
        # paper_bgcolor = colors["background"],
        xaxis=dict(
            showticklabels=True,
            tickfont={
                "size": 15
            },
            ticklen=10,
            ticks="outside",
            tickcolor=colors["background"],
            fixedrange=True
        ),
        yaxis=dict(
            range=[0, 1.2 * max(list(expression_by_tissue.values()))],
            showgrid=True,
            zeroline=True,
            # dtick=5,
            gridcolor='#f4f4f4',
            gridwidth=0.5,
            zerolinecolor='#f4f4f4',
            zerolinewidth=2,
            title="median TPM",
            title_standoff=5,
            fixedrange=True
        ),
        margin={"r": 15, "l": 35, "t": 20, "b": 10},
        legend={"traceorder": "normal",
                "x": 0.8,
                "y": 1}
    )

    return fig


# creates a combined box and bar plot from expression data

def get_bar_box_plot(expression_by_tissue):
    brain_x = [x for x in expression_by_tissue.keys() if "Brain" in x]
    brain_y = [round(expression_by_tissue[x], 3) for x in brain_x]
    other_x = [x for x in expression_by_tissue.keys() if x not in brain_x]
    other_y = [round(expression_by_tissue[x], 3) for x in other_x]

    # formatting labels
    try:
        brain_x = [re.sub(r" - ", " ", x.split("Brain - ")[1]) for x in brain_x]
        other_x = [re.sub(r" - ", " ", x) for x in other_x]
    except IndexError:  # this is due to different formats depending on source (file or gtex-api)
        brain_x = [re.sub(r"_", " ", x.split("Brain_")[1]) for x in brain_x]
        other_x = [re.sub(r"_", " ", x) for x in other_x]

    # calculating mean expression values
    brain_mean = mean(brain_y)
    other_mean = mean(other_y)

    # positioning the plots
    fig = make_subplots(rows=1, cols=3, shared_yaxes=True,
                        specs=[[{}, {"colspan": 2}, None]])

    ## all plots are being added
    fig.add_trace(go.Box(
        y=brain_y,
        text=brain_x,
        hovertemplate="<b>%{text}</b><br>%{y}",
        hoverlabel={
            "bgcolor": colors["brain"]},
        name="Brain",
        quartilemethod="linear",
        line_color=colors["brain"],
        boxpoints="all",
        showlegend=False,
        pointpos=0,
        marker={
            "color": colors["light_grey"],
            "size": 6,
            "line": {
                "color": colors["dark_grey"],
                "width": 1
            }
        },
    ),
        row=1, col=1)

    fig.add_trace(go.Box(
        y=other_y,
        text=other_x,
        hovertemplate="<b>%{text}</b><br>%{y}",
        hoverlabel={
            "bgcolor": colors["other"]},
        name="Other tissues",
        quartilemethod="linear",
        line_color=colors["other"],
        boxpoints="all",
        showlegend=False,
        pointpos=0,
        marker={
            "color": colors["light_grey"],
            "size": 6,
            "line": {
                "color": colors["dark_grey"],
                "width": 1
            }
        },
    ),
        row=1, col=1)

    fig.add_trace(go.Bar(
        x=brain_x,
        y=brain_y,
        marker={"color": colors["brain"]},
        name="Brain",
        visible=True,
        hovertemplate="<b>%{x}</b><br>%{y}",
        hoverlabel={
            "bgcolor": colors["brain"]},
    ),
        row=1, col=2)

    fig.add_trace(go.Bar(
        x=other_x,
        y=other_y,
        marker={"color": colors["other"]},
        name="Other tissues",
        visible=True,
        hovertemplate="<b>%{x}</b><br>%{y}",
        hoverlabel={
            "bgcolor": colors["other"]},
    ),
        row=1, col=2)

    # mean line brain
    fig.add_trace(go.Scatter(
        x=brain_x + other_x,
        y=[brain_mean] * len(expression_by_tissue),
        marker={"color": colors["brain"]},
        mode="lines",
        name="Brain - mean",
        visible="legendonly",
        hovertemplate="<b>Brain Mean</b>",
        hoverlabel={
            "bgcolor": colors["brain"]}
    ),
        row=1, col=2)

    # mean line other
    fig.add_trace(go.Scatter(
        x=brain_x + other_x,
        y=[other_mean] * len(expression_by_tissue),
        marker={"color": colors["other"]},
        mode="lines",
        name="Other - mean",
        visible="legendonly",
        hovertemplate="<b>Other Tissue Mean</b>",
        hoverlabel={
            "bgcolor": colors["other"]},
    ),
        row=1, col=2)

    fig.update_xaxes(showticklabels=False, row=1, col=1, fixedrange=True)
    fig.update_xaxes(showticklabels=False, row=1, col=2, fixedrange=True)

    # formats y_axis and box plot grid
    fig.update_yaxes(
        range=[0, 1.1 * max(list(expression_by_tissue.values()))],
        showgrid=True,
        zeroline=True,
        fixedrange=True,
        # dtick=5,
        gridcolor='#f4f4f4',
        gridwidth=0.5,
        zerolinecolor='#f4f4f4',
        zerolinewidth=2,
        title="median TPM",
        row=1, col=1)

    # formats bar plot grid
    fig.update_yaxes(
        range=[0, 1.1 * max(list(expression_by_tissue.values()))],
        showgrid=True,
        zeroline=True,
        fixedrange=True,
        # dtick=5,
        gridcolor='#f4f4f4',
        gridwidth=0.5,
        zerolinecolor='#f4f4f4',
        zerolinewidth=2,
        row=1, col=2)

    fig.update_layout(
        plot_bgcolor=colors["background"],
        margin={"r": 15, "l": 30, "t": 25, "b": 10},
        legend={"traceorder": "normal",
                "x": 1,
                "y": 1},
        height=400,
    )

    return fig


# parses input to achieve more tolerance to user input mistakes, checks if formatted input fits requirements
def parse_input(input):
    input = re.sub(r"^[\W]", "", input)
    input = re.sub(r"Chr|chr", "", input)
    input = re.sub(r"\s+", " ", input)
    input = re.sub(r"[^A-Za-z0-9](?!$)", ":", input)
    does_match = re.fullmatch(r"(\d{1,2}|X):\d+:[a-zA-Z]+:([a-zA-Z]*|-)+", input)
    return input, does_match


# def rate_pubtator(gene_annotation_result):


###### LAYOUT #######
app.title = "AutoCaSc"
app.layout = html.Div(
    id="app-container",
    children=[
        # ToDo implement banner
        # Banner
        # html.Div(
        #     id="banner",
        #     className="banner",
        #     children=[html.B("CAVE! This site is still under construction and is not yet fully functional!")],
        #     style={
        #         "backgroundColor": colors["red"]
        #     }
        # ),

        # changes displayed url
        dcc.Location(id='url', refresh=False),

        # stores data between callbacks
        dcc.Store(id="explanation_memory"),
        dcc.Store(id="excel_output"),
        dcc.Store(id="variant_memory"),

        # placeholder for javascript callback
        html.Div(id="copy_output"),

        # Left column
        html.Div(
            id="left-column",
            # className="four columns",
            children=[generate_description_card(), generate_control_card()],
        ),

        # Right column
        html.Div(
            id="right-column",
            children=[
                html.Div(
                    id="grid_container",
                    children=[
                        # status bar
                        html.Div(
                            id="status_bar",
                            style={
                                "backgroundColor": colors["light_grey"]
                            },
                        ),

                        # CaSc calculation results card
                        html.Div(
                            id="casc_results_card",
                            style={"min-width": "100%"},
                            hidden=True,
                            children=[
                                html.Div(
                                    html.B("Candidate Score", className="card_title"),
                                    className="card_header"
                                ),
                                html.Div(
                                    id="casc_results",
                                    children=[
                                        get_table("calculation_results_table", {"": ""}, {"": ""}, {"": ""})
                                    ]
                                ),
                                html.Div(
                                    id="copy_results_field_container"
                                )
                            ],
                        ),

                        # variant annotation results
                        html.Div(
                            html.Div(
                                id="variant_card",
                                className="pretty_container",
                                hidden=True,
                                children=[
                                    html.Div([
                                        html.B("Variant Annotation", className="card_title"),
                                    ],
                                        className="card_header"
                                    ),
                                    html.Div(id="variant_annotation",
                                             children=[
                                                 get_table("empty_variant", {"": ""}, {"": ""}, {"": ""})
                                             ]
                                             ),
                                    # html.Hr(),
                                ],
                            ),
                        ),

                        # tissue expression graph
                        html.Div(
                            id="tissue_expression_card",
                            className="pretty_container",
                            hidden=True,
                            children=[
                                html.Div([
                                    html.B("Tissue Expression", className="card_title"),
                                ],
                                    className="card_header"
                                ),
                                html.Hr(),
                                dcc.Graph(id="tissue_expression_graph", config={'displayModeBar': False}),
                            ],
                        ),
                    ],
                ),
            ]
        )
    ],
)


###### CALLBACKS ######
# changes URL to user_input
@app.callback(Output('url', 'pathname'),
              [Input("search_button", "n_clicks"),
               Input("text_input", "n_submit")],
              [State('text_input', 'value')])
def update_url(n_clicks, n_submit, input_value):
    if input_value != None:
        return "/" + input_value.replace("\n", "&")


# creates second inheritance input option depending on selected dropdown option
@app.callback(
    Output("family_history_radio_card", "children"),
    [Input("inheritance_dropdown", "value")],
    [State('text_input', 'value'),
     State('url', 'pathname')])
def update_inheritance_input(inheritance_value, text_input, url_path):
    if inheritance_value == "comphet":
        if text_input != None:
            input = text_input
        else:
            input = unquote(url_path)
        if "\n" in input or "&" in input:
            # in case of multiple lines --> impact of second allele is checked automatically
            return [
                dcc.RadioItems(
                    id="family_history_radio",
                    options=[{"label": "one affected child", "value": False},
                             {"label": "multiple affected children", "value": True}],
                    value=False),
                dcc.Dropdown(
                    id="other_impact",
                    value="auto_score",
                    style={'display': 'none'}
                ),
                html.Hr(),
                html.Hr(style={"margin-bottom": "1rem"})
            ]
        else:
            # in case of single variant
            return [
                dcc.RadioItems(
                    id="family_history_radio",
                    options=[{"label": "one affected child", "value": False},
                             {"label": "multiple affected children", "value": True}],
                    value=False),
                html.Hr(),
                html.P("Impact of second allele:", style={"fontWeight": "bold",
                                                          "margin-top": "1rem"}),
                dcc.Dropdown(
                    id="other_impact",
                    options=[
                        {"label": "high", "value": "high"},
                        {"label": "moderate", "value": "moderate"},
                        {"label": "low or modifier", "value": "low"}
                    ],
                    value="moderate",
                ),
                html.Hr(style={"margin-bottom": "1rem"})
            ]

    if inheritance_value == "homo":
        return [
            dcc.RadioItems(
                id="family_history_radio",
                options=[{"label": "one affected child", "value": False},
                         {"label": "multiple affected children", "value": True}],
                value="one_child"),
            dcc.RadioItems(
                id="other_impact"
            ),
            html.Hr(style={"margin-bottom": "1rem"})
        ]

    elif inheritance_value == "x_linked":
        return [
            dcc.RadioItems(
                id="family_history_radio",
                options=[{"label": "one boy affected", "value": False},
                         {"label": "another affected maternal male relative", "value": True}],
                value="one_boy",
                style={"margin-top": "0rem"}),
            dcc.RadioItems(
                id="other_impact"
            ),
            html.Hr(style={"margin-bottom": "1rem"})
        ]
    else:
        return [
            dcc.RadioItems(
                id="family_history_radio"
            ),
            dcc.RadioItems(
                id="other_impact"
            )
        ]


def single_variant_update(variant_of_interest, assembly, inheritance_dd_value, family_history, other_impact):
    # get variant annotation and gene scores
    variant_instance = AutoCaSc(variant_of_interest, inheritance_dd_value, family_history, other_impact, assembly)

    # status code 200 = all good, 201 = reference sequences don't match
    if variant_instance.status_code in [200, 201]:
        status_bar_color, status_bar_text, casc_card_border_color = get_status_bar(variant_instance.status_code)
    else:
        return clear_page(variant_instance.status_code, variant_of_interest)

    updated_variant_table, updated_prediction_table, casc_results_children, \
    excel_parsed, updated_graph, gtex_status_code = get_output_elements(variant_instance)

    grid_container_children = [
        # status bar
        html.Div(
            id="status_bar",
            children=status_bar_text,
            style=status_bar_color,
        ),

        # CaSc calculation results card
        html.Div(
            id="casc_results_card",
            style=casc_card_border_color,
            hidden=False,
            children=[
                html.Div(
                    html.B("Candidate Score", className="card_title"),
                    className="card_header"
                ),
                html.Div(
                    id="casc_results",
                    children=casc_results_children
                ),
                html.Div(
                    id="copy_results_field_container"
                )
            ],
        ),

        # variant annotation results
        html.Div(
            html.Div(
                id="variant_card",
                className="pretty_container",
                hidden=False,
                children=[
                    html.Div([
                        html.B("Variant Annotation", className="card_title"),
                    ],
                        className="card_header"
                    ),
                    html.Div(
                        id="variant_annotation",
                        children=[updated_variant_table]
                    ),
                    # html.Hr(),
                ],
            ),
        ),

        # tissue expression graph
        html.Div(
            id="tissue_expression_card",
            className="pretty_container",
            hidden=False,
            children=[
                html.Div([
                    html.B("Tissue Expression", className="card_title"),
                ],
                    className="card_header"
                ),
                html.Hr(),
                dcc.Graph(id="tissue_expression_graph",
                          config={'displayModeBar': False},
                          figure=updated_graph),
            ],
        ),
    ]

    return ["",  # empty output for loading sign
            variant_of_interest,
            variant_instance.explanation_dict,
            excel_parsed,
            grid_container_children
            ]


# def batch_update(variant_input, assembly, inheritance_dd_value, family_history):
#     variant_list = variant_input.split("&")
#     variant_dict = {}
#     any_success = False
#
#     for variant in variant_list:
#         autocasc_instance = AutoCaSc(variant, inheritance_dd_value, family_history, assembly=assembly)
#         if autocasc_instance.status_code == 200:
#             variant_dict[variant] = autocasc_instance
#             any_success = True
#         else:
#             variant_dict[variant] = None
#
#     if any_success == False:
#         return clear_page(400)
#     else:
#         status_bar_color, status_bar_text, casc_card_border_color = get_status_bar(200)
#
#     # this module rescores variants taking the impact of the second allele into account
#     if inheritance_dd_value == "comphet":
#         gene_dict = {}
#         for variant in variant_dict.keys():
#             gene_dict[variant] = variant_dict.get(variant).gene_id
#         comphet_df = pd.DataFrame.from_dict(gene_dict, orient='index').reset_index()
#         comphet_df.columns = ['variant', 'gene_id']
#
#         for iter_var in comphet_df.variant:
#             try:
#                 variant_gene = comphet_df.loc[comphet_df.variant == iter_var, "gene_id"].values[0]
#                 other_variant = comphet_df.loc[comphet_df.gene_id == variant_gene, "variant"].loc[comphet_df.variant != iter_var].values[0]
#                 other_impact = variant_dict.get(other_variant).impact
#                 variant_dict[iter_var].other_impact = other_impact
#                 variant_dict[iter_var].get_scores()
#
#             except IndexError:
#                 pass
#
#     clipboard_list = []
#     for variant in variant_list:
#         clipboard_list.append(get_clipboard_format_single(variant_dict[variant].__dict__))
#     clipboard_data = {"data": "\n".join(clipboard_list)}
#
#     casc_results_children = get_casc_batch_table(variant_dict)
#
#     grid_container_children = [
#                         # status bar
#                         html.Div(
#                             id="status_bar",
#                             children=status_bar_text,
#                             style=status_bar_color,
#                         ),
#
#                         # CaSc calculation results card
#                         html.Div(
#                             id="casc_results_card",
#                             style=casc_card_border_color,
#                             hidden=False,
#                             children=[
#                                 html.Div(
#                                     html.B("Candidate Score", className="card_title"),
#                                     className="card_header"
#                                 ),
#                                 html.Div(
#                                     id="casc_results",
#                                     children=casc_results_children
#                                 ),
#                                 html.Br(),
#                                 html.Div(
#                                     [
#                                         html.Button(["copy results"],
#                                                     id="copy_button",
#                                                     className="six columns"),
#                                     ]
#                                 ),
#                             ],
#                         ),
#             ]
#
#     return ["",  # empty output for loading sign
#             "\n".join(variant_list),
#             {},  # ToDo explanation dicts
#             clipboard_data,  # ToDo excel_parsed
#             grid_container_children
#             ]


# core program, calling aquisition of data, calculation of CaSc and returns output data
def batch_update(variant_input, assembly, inheritance_dd_value, family_history):
    variant_list = variant_input.split("&")
    variant_dict = {}
    any_success = False

    for variant in variant_list:
        autocasc_instance = AutoCaSc(variant, inheritance_dd_value, family_history, assembly=assembly)
        if autocasc_instance.status_code == 200:
            variant_dict[variant] = autocasc_instance
            any_success = True
        else:
            variant_dict[variant] = None

    if any_success == False:
        return clear_page(400)
    else:
        status_bar_color, status_bar_text, casc_card_border_color = get_status_bar(200)

    # this module rescores variants taking the impact of the second allele into account
    if inheritance_dd_value == "comphet":
        gene_dict = {}
        for variant in variant_dict.keys():
            gene_dict[variant] = variant_dict.get(variant).gene_id
        for _variant in gene_dict.keys():
            try:
                variant_gene = gene_dict.get(_variant)
                other_variant = list(filterTheDict(gene_dict, variant_gene, _variant).keys())[0]
                other_impact = variant_dict.get(other_variant).impact
                variant_dict[_variant].other_impact = other_impact
                variant_dict[_variant].get_scores()
            except IndexError:
                pass

        comphet_id = 0  # unique identifier for clear identification of corresponding compound heterozygous variants
        for _variant in gene_dict.keys():
            if variant_dict.get(_variant).comphet_id:
                continue
            variant_gene = gene_dict.get(_variant)
            other_variant = list(filterTheDict(gene_dict, variant_gene, _variant).keys())[0]

            combined_variant_score = mean([variant_dict.get(_variant).variant_score,
                                           variant_dict.get(other_variant).variant_score])
            variant_dict.get(_variant).variant_score = combined_variant_score
            variant_dict.get(other_variant).variant_score = combined_variant_score
            variant_dict.get(_variant).comphet_id = comphet_id
            variant_dict.get(other_variant).comphet_id = comphet_id
            variant_dict.get(_variant).calculate_canidate_score()
            variant_dict.get(other_variant).calculate_canidate_score()
            comphet_id += 1

    clipboard_list = []
    for variant in variant_list:
        clipboard_list.append(get_clipboard_format_single(variant_dict[variant].__dict__))
    clipboard_data = {"data": "\n".join(clipboard_list)}

    casc_results_children = get_casc_batch_table(variant_dict)

    grid_container_children = [
        # status bar
        html.Div(
            id="status_bar",
            children=status_bar_text,
            style=status_bar_color,
        ),

        # CaSc calculation results card
        html.Div(
            id="casc_results_card",
            style=casc_card_border_color,
            hidden=False,
            children=[
                html.Div(
                    html.B("Candidate Score", className="card_title"),
                    className="card_header"
                ),
                html.Div(
                    id="casc_results",
                    children=casc_results_children
                ),
                html.Br(),
                html.Div(
                    [
                        html.Button(["copy results"],
                                    id="copy_button",
                                    className="six columns"),
                    ]
                ),
            ],
        ),
    ]

    return ["",  # empty output for loading sign
            "\n".join(variant_list),
            {},  # ToDo explanation dicts
            clipboard_data,  # ToDo excel_parsed
            grid_container_children
            ]


@app.callback(
    # list cointaining all the elements that are being updated by the function by returning corresponding values
    [Output("loading-output-1", "children"),
     Output("text_input", "value"),
     Output("explanation_memory", "data"),
     Output("excel_output", "data"),
     Output("grid_container", "children")
     ],

    # trigger for updated function
    [Input("url", "pathname")],

    # values that are needed for updating but don't trigger the update
    [State("dd_reference_genome", "value"),
     State("inheritance_dropdown", "value"),
     State("family_history_radio", "value"),
     State("other_impact", "value")])
def update(variant_input, assembly, inheritance_dd_value, family_history, other_impact):
    """
    The core function of the web app coordinating data updating, formatting and error handling.
    :param variant_input: variant genomic position, reference and alternative sequences
    :param assembly: GRCh37 or GRCh38
    :param inheritance_dd_value: selected inheritance mode of variant
    :param family_history: are there multiple affected family members? either True or False
    :param other_impact: in case of compund heterozygous variants, impact of other variant
    :return:
    """
    if variant_input != "":
        variant_input = unquote(variant_input).strip("/&")
    # ToDo variant_input ändern in variant_pos oder so
    if variant_input != "":
        if "&" in variant_input:
            return batch_update(variant_input, assembly, inheritance_dd_value, family_history)
        else:
            return single_variant_update(variant_input, assembly, inheritance_dd_value, family_history, other_impact)

    else:
        return clear_page(300)


app.clientside_callback(
    """
    function placeholder(n_clicks, excel_output) {
        const copyToClipboard = str => {
            const el = document.createElement('textarea');
            el.value = str;
            el.setAttribute('readonly', '');
            el.style.position = 'absolute';
            el.style.left = '-9999px';
            document.body.appendChild(el);
            const selected =
            document.getSelection().rangeCount > 0 ? document.getSelection().getRangeAt(0) : false;
            el.select();
            document.execCommand('copy');
            document.body.removeChild(el);
            if (selected) {
            document.getSelection().removeAllRanges();
            document.getSelection().addRange(selected);
            }
        };

        var userLang = navigator.language || navigator.userLanguage;
        var output_text = excel_output.data;

        if (userLang.includes("de")) {
            output_text = output_text.replace(/\./g, ",")
        }

        copyToClipboard(output_text);

        }
    """,
    [Output("copy_button", "children")],
    [Input("copy_button", "n_clicks")],
    [State("excel_output", "data")]
)

# @app.callback(
#     Output("copy_button", "children"),
#     [Input("copy_button", "n_clicks")],
#     [State("copy_button_memory", "data")])
# def button_updated(n_clicks, data):
#     data["clicks"] = data.get("clicks") + 1
#     return data.get("clicks")
#
# @app.callback(
#     Output("copy_button_memory", "data"),
#     [Input("copy_button", "children")]
# )
# def temp_function(n_clicks):
#     return {"clicks": n_clicks}

##### RUN #####
if __name__ == '__main__':
    app.run_server(debug=False, dev_tools_hot_reload=True)
