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
from numpy import array

from AutoCaSc_core.AutoCaSc import AutoCaSc
from dash_extensions import Download

server = Flask(__name__)
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP],
                server=server)
app.title = "webAutoCaSc"

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

url_bar_and_content_div = html.Div([
    dcc.Location(id='url', refresh=False),
    dcc.Store(id="query_memory"),
    dcc.Store(id="variant_queue_input"),
    dcc.Store(id="variant_queue_url"),
    dcc.Store(id="variant_memory"),
    dcc.Store(id="results_memory"),
    Download(id="download"),
    navbar,
    html.Div(id='page-content'),
    footer
])

variant_input_card = dbc.FormGroup(
    [
        dbc.Input(
            type="text",
            id="variant_input",
            placeholder="e.g. X:12345:T:C",
            autoFocus=True
        ),
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
    html.Div(style={"height": "10vh"}),
    dbc.Container(
            [
                dbc.Row(
                    dbc.Col(
                        [
                            dcc.Markdown("""# Welcome to **webAutoCaSc**,\n
#### a webinterface for the automatic CaSc classification of research candidate variants in neurodevelopmental disorders."""),
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
            style={"max-height": "calc(100vh - 150px)",
                    "overflow-y": "auto"}
        ),
    # html.Div(style={"height": "1vh"}),
])


search_page = html.Div([
    dbc.Container([
        dbc.Spinner(fullscreen=True)
    ])
])

results_page_clear = html.Div([
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
        style={"max-height": "calc(100vh - 150px)",
               "overflow-y": "auto"}
    ),
],

)
citations = html.Div([
    html.Hr(),
    dbc.Container([
        html.P([
            "1. ", html.A(html.B("VEP"), href="https://grch37.ensembl.org/info/docs/tools/vep/index.html", target="_blank"),
            ": McLaren, W. et al. The Ensembl Variant Effect Predictor. Genome Biol 17, 122 (2016)."
        ]),
        html.P([
            "2. ", html.A(html.B("gnomAD"), href="https://gnomad.broadinstitute.org", target="_blank"),
            ": Karczewski, K. J. et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 581, 434–443 (2020)."
        ]),
        html.P([
            "3. ", html.A(html.B("GTEx"), href="http://www.gtexportal.org/home/index.html", target="_blank"),
            ": Consortium, T. Gte. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369, 1318–1330 (2020)."
        ]),
        html.P([
            "4. ", html.A(html.B("STRING"), href="https://string-db.org", target="_blank"),
            ": Szklarczyk, D. et al. STRING v11: protein-protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets. Nucleic Acids Res 47, D607–D613 (2019)."
        ]),
        html.P([
            "5. ", html.A(html.B("MGI"), href="http://www.informatics.jax.org", target="_blank"),
            ": Bult, C. J. et al. Mouse Genome Database (MGD) 2019. Nucleic Acids Res 47, D801–D806 (2019)."
        ]),
        html.P([
            "6. ", html.A(html.B("PubTator"), href="https://www.ncbi.nlm.nih.gov/research/pubtator/", target="_blank"),
            ": Wei, C.-H., Allot, A., Leaman, R. & Lu, Z. PubTator central: automated concept annotation for biomedical full text articles. Nucleic Acids Res 47, W587–W593 (2019)."
        ]),
        html.P([
            "7. ", html.A(html.B("PsyMuKB"), href="http://www.psymukb.net", target="_blank"),
            ": Lin, G. N. et al. PsyMuKB: An Integrative De Novo Variant Knowledge Base for Developmental Disorders. Genomics Proteomics Bioinformatics 17, 453–464 (2019)."
        ]),
        html.P([
            "8. ", html.A(html.B("DisGeNET"), href="https://www.disgenet.org", target="_blank"),
            ": Piñero, J. et al. DisGeNET: a comprehensive platform integrating information on human disease-associated genes and variants. Nucleic Acids Res 45, D833–D839 (2017)."
        ])
    ]),
    html.Hr()
])

about_eng = [html.P("AutoCaSc is a tool for quantifying the plausibility of candidate variants for Neurodevelopmental Disorders (NDD). AutoCaSc is intended to be used on impactful rare variants in a research setting. In its current version, 12 parameters are counted in, achieving a maximum of 15 points. User inputs are the identified variant in a standard HGVS/VCF format together with segregation aspects (de novo, recessive, dominant and X-chromosomal). We use the Ensembl VEP REST API (1) to annotate variant attributes (e.g. variant consequence, allele frequency from gnomAD, in silico predictions) and gene based scores dependent on inheritance mode (e.g. high Z-score is of relevant for de novo missense) from dbNSFP. Other attributes were previously labor intensive and predisposed to variability. These included important categories like expression in the nervous system, neuronal functions, co-expression and protein interactions, search for relevant literature, model organisms and observations in screening studies. As an objective approach we now searched a defined set of databases (gnomAD (2), GTEx (3), STRING (4), MGI (5), PubTator (6), PsyMuKB (7), DisGeNET (8)) and generated empirical cut-offs for each category by comparing the respective readout between a manually curated list of known NDD genes from the SysID database (9) and a list of genes not involved in NDD.",
                   style={"text-align": "justify"}),
            html.Br(),
            html.P("Feel free to contact johann.lieberwirth@medizin.uni-leipzig.de or rami.aboujamra@medizin.uni-leipzig.de in case you have further questions or in case you have found a bug.",
                   style={"text-align": "justify"})]

about_ger = [html.P("AutoCaSc ist ein Skript zum automatisierten Bewerten von Kandidatenvarianten in Fällen neuronaler Entwicklungsverzögerung. Es ist ausschließlich für Forschungszwecke zu benutzen. Die Annotation der Varianten erfolgt mit der REST API von VEP (ensembl, (1)). Zur Berechnung der Kandidatenpunktzahl (Candidate score, CaSc) werden 12 verschiedene Parameter einbezogen. Dies sind die Art der Vererbung (z.B. de novo), Genattribute wie pLI & Z-Score (gnomAD (2)), Expressionsmuster (GTEx (3)), in silico Analysen, Proteininteratkionsdatnebanken (StringDB (4)), Tierdatenbanken (MGI (5)), Literaturdatenbanken (Pubtator Central (6)),  sowie weitere (PsymuKB (7), DisGeNET (8)). Die maximal erreichbare Punktzahl sind 15 Punkte. Je höher der erreichte Punktwert, desto plausibler scheint die aus den zugrundeliegenden Daten errechnete Pathogenität der Variante mit Blick auf neuronale Entwicklungsverzögerung.",
                   style={"text-align": "justify"}),
             html.Br(),
             html.P("Bei Fragen und Anmerkungen kontaktieren Sie bitte johann.lieberwirth@medizin.uni-leipzig.de oder rami.aboujamra@medizin.uni-leipzig.de.",
                   style={"text-align": "justify"})]

about_page = html.Div([
    html.Br(),
    dbc.Container(
        [
            dbc.Row([
                dbc.Col(html.H2("About"),
                        width="auto"),
                dbc.Col(dbc.Button("DE", id="about_language_button"),
                        width="auto")
            ]),
            html.Br(),
            html.Div(about_eng,
            id="about_text"),
            html.Br(),
            citations
        ],
        style={"max-height": "calc(100vh - 180px)",  # html.Hr margin seemed to induce global scrollbar
               "overflow-y": "auto"}
    ),
    # html.Div(style={"height": "1vh"}),
    ])


faq_ger = dcc.Markdown("""
            __Was ist AutoCaSc?__  
            AutoCaSc ist ein Werkzeug zur systematischen Evaluierung der Plausibilität von Varianten in Genen, welche bislang nicht mit Erkrankungen in Verbindung gebracht wurden ("Kandidatengene"), in Fällen neurologiscer Entwicklungsverzögerung (neurodevelopmental disorder, NDD). Solche Varianten werden üblicherweise durch genomweites Screeningmethoden bei Individuen mit NDD, aber ohne eindeutige diagnostische Variante identifiziert. AutoCaSc berücksichtigt variantenspezifische Parameter (Konservierung, "in silico"-Vorhersagen), genspezifische Parameter ("gene constraint", Expressionsmuster, Proteininteraktionen), die Segregation der Variante und das Gesamtzusammenspiel zwischen diesen Parametern.

            __Wofür stehen die 4 Unterscores?__  
            - __Variant Attributes (6 Punkte max):__ Dazu gehören Konservierung (GERP++), "in silico"-Vorhersagen (MutationAssessor, MutationTaster, Sift), Spleißstellenvorhersagen (MaxEntScan, AdaBoost, RandomForest) und erwartete Auswirkungen (VEP).
            - __Gene Attributes (1 Punkte max):__ Es handelt sich dabei um gene-constraint Parameter aus gnomAD; LOUEF für Loss-of-Function-Varianten, Z für Missense-Varianten.
            - __Inheritance (2 Punkte max):__ Diese Punkte hängen von der Vererbung und der Segregation der Variante in der Familie ab.
            - __Gene Plausibility (6 Punkte max):__ Diese Punkte werden auf der Grundlage des Expressionsmusters des Gens, der Protein-Protein-Interaktionen, der Phänotypen in Tiermodellen, der in PubMed veröffentlichten Artikel zum Gen, de novo-Varianten im Gen die mit NDD in Verbindung gebacht wurden, und anderer Quellen berechnet.

            __Wie kann man mehrere (compound heterozygote) Variaten eingeben?__  
            Mehrere Varianten können eingegeben werden, indem sie durch ein Komma getrennt werden. Wenn "compound heterozygous" ausgewählt ist, findet webAutoCaSc automatisch Varianten im gleichen Gen und verarbeitet diese als entsprechende compound heterozygote Varianten.

            __Wofür stehen die Vererbungsoptionen?__  
            - __De novo:__ De novo-Varianten werden nur im Index identifiziert und sind nicht von den Eltern vererbt worden.
            - __Inherited dominant:__ Im Falle einer vererbten dominanten Variante wurde die Variante von einem ebenfalls betroffenen Elternteil vererbt.
            - __Homozygous recessive:__ Die Variante wird homozygot im Index und in heterozygot in beiden gesunden Elternteilen identifiziert.
            - __X-linked:__ X-chromosomale Varianten werden von der heterozygoten Trägermutter vererbt und verursachen bei einem männlichen Nachkommen einen Phänotyp, da er nur ein betroffenes Allel und kein gesundes zweites Allel zum Ausgleich hat. De-novo-Varianten auf dem X-Chromosom, sowohl bei weiblichen als auch bei männlichen Index-Individuen, werden in der Option De-novo-Vererbung berücksichtigt.
            - __Compound heterozygous:__ Compound heterozygote Varianten sind zwei verschiedene Varianten im selben Gen, aber auf verschiedenen Allelen. Jede wird von nur einem heterozygoten Trägerelternteil vererbt.
            - __Unknown:__ Die Option "unbekannt" kann verwendet werden, wenn Informationen zu den Eltern und damit zur Segregation fehlen.
            
            __Ab wann ist ein Score hoch?__
            Die maximale Punktzahl sind 15. Um ein besseres Gefühl zu geben, ob ein candidate score hoch ist, wurde ein Tooltip implementiert, welcher angezeigt wird wenn der Mauszeiger über dem Ergebnis schwebt. Der Tooltip zeigt an, wie viel Prozent der Kandidaten, welche am Institut für Humangenetik in Leipzig evaluiert wurden, ein höheres Ergebnis erreichten.

            __Wofür steht _webAutoCaSc_?__  
            _AutoCaSc_ steht für __Auto__mated __Ca__ndidate __Sc__ore. Wir benutzen den __web__ Präfix, um das command line interface (CLI) von der Webapp zu unterscheiden, welche auf dem AutoCaSc script basiert. Das __Ca__ndidate __Sc__ore Prinzip wurde bereits von _Büttner et al. bioRxiv. 2019_ beschrieben. 

            __Kann webAutoCaSc für andere Phänotypen genutzt werden?__  
            AutoCaSc wurde für die Arbeit mit NDDs entwickelt. Wir empfehlen nicht, es für andere Phänotypen zu verwenden. Für zukünftige Updates planen wir ein verallgemeinertes, phänotyp-unabhängiges Framework.
            """)

faq_eng = dcc.Markdown("""
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
            
            __When is a score high?__  
            The maximum score is 15. To give a better feeling if a candidate score is high, a tooltip has been implemented, which is displayed when the mouse pointer hovers over the result. The tooltip shows what percentage of candidates evaluated at the Institute of Human Genetics in Leipzig achieved a higher score.
            
            __What does _webAutoCaSc_ stand for?__  
            _AutoCaSc_ stands for __Auto__mated __Ca__ndidate __Sc__ore. We use the __web__ prefix to distinguish the command line interface (CLI) from the webapplication running the AutoCaSc script. The __Ca__ndidate __Sc__ore principle has been previously described (Büttner et al. bioRxiv. 2019). 
            
            __Can webAutoCaSc be used for other phenotypes as well?__  
            AutoCaSc has been developed to work for NDDs. We don't recommend using it for other phenotypes. We are planning a generalized phenotype agnostic framework for future updates.
            """)


faq_page = html.Div([
    dbc.Container(
        [
            html.Br(),
            dbc.Row([
                dbc.Col(html.H2("FAQ"),
                        width="auto"),
                dbc.Col(dbc.Button("DE", id="faq_language_button"),
                        width="auto")
            ]),
            html.Br(),
            html.Div(faq_eng,
                     id="faq_text")
        ],
        style={"max-height": "calc(100vh - 150px)",
               "overflow-y": "auto"}
    ),
    ])

impressum_ger = dcc.Markdown("""
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

impressum_eng = [
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
            ]


impressum_page = html.Div(
    [
        dbc.Container([
            html.Br(),
            dbc.Row([
                dbc.Col(html.H2("Impressum"),
                        width="auto"),
                dbc.Col(dbc.Button("EN", id="impressum_language_button"),
                        width="auto")
            ]),
            html.Br(),
            html.Div(impressum_ger,
            id="impressum_text",
            )
        ],
        style={"max-height": "calc(100vh - 150px)",
               "overflow-y": "auto"}
        ),
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
    faq_page
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
                style={"max-height": "calc(100vh - 150px)",
                       "overflow-y": "auto"}
        ),
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


def get_percentile(candidate_score):
    try:
        int(candidate_score)
        internal_scores = [10.46, 9.24, 9.86, 9.51, 11.81, 9.24, 10.5, 8.42, 8.03, 9.56, 0.0, 10.52, 8.04, 9.29, 9.66, 9.02, 7.22, 8.33, 10.37, 9.97, 0.0, 9.06, 6.87, 10.2, 9.06, 9.35, 8.01, 8.63, 8.94, 7.58, 8.47, 8.12, 8.9, 8.45, 9.19, 8.34, 5.85, 5.02, 7.96, 6.69, 8.14, 0.0, 7.63, 7.81, 8.06, 7.62, 6.28, 8.81, 8.41, 10.34, 9.64, 0.0, 8.23, 5.99, 0.0, 8.09, 8.78, 9.57, 6.61, 0.0, 8.43, 9.35, 5.74, 6.84, 6.39, 7.54, 9.04, 9.23, 9.31, 6.35, 10.0, 4.82, 7.9, 5.08, 8.5, 6.6, 8.05, 5.95, 6.2, 8.77, 6.49, 7.53, 0.0, 6.75, 0.0, 6.46, 5.0, 8.05, 7.76, 5.5, 5.79, 4.36, 8.76, 6.31, 5.86, 5.07, 7.53, 0.0, 6.79, 0.0, 9.53, 8.09, 6.63, 6.11, 7.24, 10.86, 5.57, 5.55, 6.73, 7.9, 7.34, 6.81, 7.66, 5.56, 6.1, 7.47, 6.03, 8.32, 3.85, 7.21, 5.46, 0.0, 8.32, 7.79, 9.16, 5.83, 5.05, 5.18, 5.41, 7.51, 6.28, 5.71, 5.99, 5.12, 6.41, 0.0, 5.88, 5.31, 6.46, 4.95, 0.0, 6.04, 5.12, 6.44, 4.74, 3.95, 4.37, 6.21, 4.66, 5.52, 5.76, 4.66, 6.24, 8.1, 5.15, 4.89, 4.94, 8.04, 0.0, 7.23, 5.02, 6.66, 4.94, 5.04, 7.97, 4.71, 5.84, 6.77, 5.67, 6.04, 4.6, 6.24, 5.85, 6.68, 5.07, 5.55, 7.57, 5.98, 7.99, 4.21, 6.89, 8.4, 6.89, 5.55, 5.18, 5.93, 6.88, 6.83, 6.22, 7.68, 5.95, 6.62, 4.05, 3.7, 6.09, 8.1, 0.0, 6.42, 4.19, 5.26, 0.0, 0.0, 4.03, 5.23, 3.96, 6.4, 5.53, 8.32, 5.92, 8.25, 0.0, 0.0, 5.7, 5.36, 4.99, 0.0, 0.0, 5.04, 6.28, 4.96, 0.0, 0.0, 5.94, 4.73, 6.16, 4.85, 7.14, 0.0, 5.05, 5.44, 0.0, 5.12, 5.84, 5.0, 6.46, 4.45, 6.75, 4.7, 0.0, 0.0, 3.45, 5.44, 6.22, 5.24, 0.0, 4.94, 5.83, 0.0, 5.13, 6.45, 4.65, 5.64, 6.02, 6.34, 5.94, 5.33, 5.04, 4.89, 6.21, 4.36, 0.0, 0.0, 6.23, 0.0, 6.22, 5.12, 5.94, 4.18, 6.03, 5.18, 5.51, 5.04, 5.49, 6.15, 3.83, 7.9, 4.88, 8.07, 6.67, 5.28, 5.26, 5.19, 4.95, 5.0, 4.26, 7.82, 5.8, 5.89, 5.98, 5.02, 5.39, 6.69, 0.0, 0.0, 6.98, 8.03, 4.63, 5.75, 6.43, 6.5, 6.6, 6.22, 4.41, 6.32, 4.57, 5.41, 5.73, 0.0, 4.82, 5.53, 0.0, 5.26, 0.0, 6.77, 4.72, 4.46, 4.66, 5.24, 4.31, 0.0, 4.6, 6.11, 4.74, 4.23, 4.75, 4.91, 4.96, 4.73, 4.58, 3.41, 4.92, 0.0, 5.45, 7.24, 4.08, 5.79, 0.0, 4.76, 5.1, 6.11, 0.0, 4.27, 4.98, 0.0, 4.17, 5.05, 4.88, 6.78, 4.73, 5.03, 5.6, 5.94, 0.0, 6.22, 4.38, 5.76, 4.35, 5.15, 5.91, 4.54, 0.0, 5.17, 0.0, 4.84, 4.8, 4.23, 5.09, 4.38, 4.62, 4.65, 4.5, 6.74, 5.51, 4.26, 5.19, 4.88, 6.01, 4.38, 5.63, 4.84, 0.0, 0.0, 5.1, 5.7, 4.98, 4.84, 6.35, 5.54, 4.09, 5.0, 0.0, 0.0, 4.01, 5.54, 3.91, 0.0, 4.64, 8.49, 3.56, 6.4, 0.0, 4.59, 4.15, 7.2, 4.79, 4.85, 0.0, 0.0, 3.7, 3.84, 4.8, 4.04, 5.14, 6.07, 5.81, 4.81, 0.0, 4.07, 4.73, 4.78, 4.32, 5.2, 4.93, 0.0, 2.28, 3.25, 3.5, 4.28, 4.78, 4.15, 4.29, 6.32, 7.4, 5.54, 4.19, 0.0, 4.84, 0.0, 5.21, 5.36, 3.71, 0.0, 4.92, 3.68, 6.24, 4.39, 5.22, 0.0, 4.3, 3.71, 3.94, 2.86, 5.36, 2.53, 0.0, 7.04, 0.0, 4.55, 5.1, 4.16, 4.53, 4.91, 5.29, 6.18, 0.0, 4.61, 6.32, 3.57, 4.85, 0.0, 3.2, 2.29, 0.0, 0.0, 0.0, 4.39, 4.01, 5.66, 5.46, 0.0, 4.37, 4.96, 3.67, 5.89, 5.01, 2.46, 3.56, 0.0, 3.83, 4.61, 5.31, 0.0, 4.58, 3.5, 6.04, 5.21, 0.0, 1.26, 2.17, 0.0, 0.0, 3.16, 4.38, 3.08, 2.33, 4.61, 4.57, 4.84, 3.14, 0.0, 3.32, 4.03, 4.8, 4.92, 4.73, 5.81, 3.12, 0.0, 5.58, 6.04, 3.55, 3.59, 0.0, 0.0, 3.77, 5.69, 3.68, 0.0, 3.4, 3.13, 0.0, 3.99, 5.45, 4.66, 4.67, 3.68, 3.77, 0.0, 0.0, 4.71, 0.0, 4.01, 3.48, 2.39, 3.73, 3.29, 5.18, 0.0, 0.0, 0.0, 3.32, 2.86, 0.0, 0.0, 5.02, 4.25, 6.38, 0.0, 0.0, 2.83, 0.0, 4.19, 0.0, 3.28, 3.45, 4.72, 3.54, 4.96, 3.12, 0.0, 4.54, 4.61, 3.44, 3.08, 1.93, 4.63, 3.1, 0.0, 3.96, 4.59, 3.26, 5.0, 5.17, 4.85, 0.0, 2.89, 3.89, 4.92, 3.32, 3.37, 0.0, 3.38, 5.64, 2.66, 3.5, 4.78, 0.0, 0.0, 4.1, 3.42, 0.0, 2.74, 0.0, 3.55, 2.08, 5.14, 0.0, 1.04, 0.0, 3.68, 2.72, 2.35, 1.69, 3.37, 0.0, 9.13, 6.9, 6.69, 2.8, 4.33, 7.13, 9.79, 5.45, 5.02, 7.67, 4.32, 4.3, 9.23, 4.52, 5.54, 7.13, 5.8, 4.15, 3.93, 4.35, 4.43, 2.95, 1.62, 4.24, 5.29, 4.82, 5.25, 5.0, 5.7, 3.18, 5.17, 3.36, 8.17, 4.94, 5.9, 4.24, 5.77, 7.12, 4.66, 2.73, 5.03, 0.0, 4.5, 4.75, 3.31, 0.0, 4.27, 3.4, 5.37, 2.59, 4.14, 4.85, 4.66, 7.03, 5.0, 5.96, 5.27, 5.44, 5.36, 6.1]

        a = array(internal_scores)
        percentile = 100. * len(a[a > candidate_score]) / len(internal_scores)
        return str(round(percentile)) + "%"
    except (TypeError, ValueError):
        return "-"


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
                        html.Th(_instance_attributes.get("candidate_score"),
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
                                    dbc.Col(html.H3(f"Candidate Score: {_instance_attributes.get('candidate_score')}",
                                                    id="percentile_target")),
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
                    dbc.Col(html.H3(f"Candidate Score: {_instance_attributes.get('candidate_score')}",
                                    id="percentile_target"),
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
                dbc.Tooltip(f"About {get_percentile(_instance_attributes.get('candidate_score'))} of the scored variants @Institute for Human Genetics Leipzig had a higher score than this.",
                            target="percentile_target"),
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
        return about_ger, "EN"
    else:
        return about_eng, "DE"


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
        return impressum_ger, "EN"
    else:
        return impressum_eng, "DE"


@app.callback(
    Output("faq_text", "children"),
    Output("faq_language_button", "children"),
    Input("faq_language_button", "n_clicks"),
    State("faq_language_button", "children")
)
def get_faq_text(n_clicks, language):
    if not n_clicks:
        raise PreventUpdate
    if language == "DE":
        return faq_ger, "EN"
    else:
        return faq_eng, "DE"

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
        df.loc[i, "candidate_score"] = _instance.get("candidate_score")
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