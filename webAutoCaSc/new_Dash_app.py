import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

url_bar_and_content_div = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
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
            # todo increase text size
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

results_page = html.Div([
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
    ])
])

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
                        color="primary"
                    )
                )
            )
        ]
    )
])

# index layout
app.layout = url_bar_and_content_div

# "complete" layout
app.validation_layout = html.Div([
    url_bar_and_content_div,
    about_page,
    landing_page,
    results_page
])


# Index callbacks
@app.callback(Output('page-content', 'children'),
              Input('url', 'pathname'))
def display_page(pathname):
    if pathname == "/about":
        return about_page
    if pathname == "/results":
        return results_page
    else:
        return landing_page



@app.callback(
        Output(f"navbar-collapse", "is_open"),
        [Input(f"navbar-toggler", "n_clicks")],
        [State(f"navbar-collapse", "is_open")])
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


if __name__ == '__main__':
    app.run_server(debug=False,
                   dev_tools_hot_reload=True)