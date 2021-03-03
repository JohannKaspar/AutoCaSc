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
        internal_scores = [13.7, 13.5, 13.33, 13.0, 12.66, 12.66, 12.2, 11.9, 11.9, 11.83, 11.7, 11.690000000000001, 11.61,
                           11.56, 11.56, 11.530000000000001, 11.5, 11.5, 11.4, 11.33, 11.2, 11.2, 10.99, 10.85,
                           10.798333333333332, 10.7, 10.7, 10.690000000000001, 10.690000000000001, 10.650323333333333,
                           10.65, 10.49, 10.4, 10.33, 10.24, 10.21, 10.2, 10.16, 10.030000000000001, 10.0,
                           9.976666666666667, 9.859743333333334, 9.83, 9.796666666666665, 9.26798, 9.763333333333334, 9.76,
                           9.73, 9.716666666666667, 9.704629999999998, 9.7, 9.7, 9.690000000000001, 9.66, 9.66, 9.605,
                           9.556666666666667, 9.55, 9.546429999999999, 9.513333333333334, 9.4, 9.4, 9.4, 9.4, 9.36, 9.36,
                           9.36, 9.33, 9.326666666666666, 9.324644999999999, 9.3, 9.280000000000001, 9.275,
                           9.230243333333332, 9.200000000000001, 9.2, 9.2, 9.2, 9.2, 9.18, 9.1165, 9.103333333333333,
                           9.06017, 9.030000000000001, 9.030000000000001, 9.01, 9.009893333333334, 9.000000000000002,
                           8.98489, 8.966666666666667, 8.946666666666667, 8.933333333333334, 8.914876666666666,
                           8.900406666666667, 8.9, 8.9, 8.873333333333333, 8.842799999999999, 8.836666666666666, 8.83, 8.81,
                           8.799999999999999, 8.740000000000002, 8.73, 8.71, 8.7, 8.7, 8.7, 8.7, 8.7, 8.7, 8.7,
                           8.693000000000001, 8.686666666666667, 8.673333333333334, 8.67, 8.66, 8.64, 8.619250000000001,
                           8.61, 8.6, 8.58, 8.57, 8.566666666666666, 8.56, 8.559693333333334, 8.530000000000001,
                           8.530000000000001, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.440000000000001,
                           8.416666666666666, 8.4, 8.4, 8.39, 8.383333333333335, 8.360000000000001, 8.357876666666666,
                           8.35585, 8.353333333333333, 8.35, 8.299999999999999, 8.287923333333334, 8.280000000000001,
                           8.276666666666667, 8.25536, 8.23, 8.23, 8.21914, 8.213043333333333, 8.206800000000001, 8.2,
                           8.173783333333333, 8.161666666666667, 8.15, 8.145, 8.131383333333334, 8.10576, 8.100000000000001,
                           8.076666666666668, 8.060783333333333, 8.06, 8.059999999999999, 8.055226666666666, 8.05,
                           8.043333333333333, 8.039976666666666, 8.030000000000001, 8.030000000000001, 8.030000000000001,
                           8.030000000000001, 8.029896666666666, 8.016666666666667, 8.013333333333334, 8.01, 7.99, 7.99,
                           7.97232, 7.970000000000001, 7.9657, 7.964833333333334, 7.933333333333334, 7.930000000000001,
                           7.906666666666666, 7.9, 7.9, 7.9, 7.9, 7.876666666666667, 7.83, 7.81925, 7.816666666666666,
                           7.796666666666667, 7.79, 7.75, 7.732243333333333, 7.73, 7.73, 7.729333333333333, 7.7, 7.7, 7.7,
                           7.7, 7.7, 7.7, 7.7, 7.686986666666667, 7.680000000000001, 7.676666666666667, 7.666666666666668,
                           7.61531, 7.6, 7.57, 7.562666666666667, 7.555893333333334, 7.55276, 7.5505700000000004,
                           7.536666666666667, 7.5166666666666675, 7.51, 7.5, 7.486666666666667, 7.483333333333333, 7.47923,
                           7.474753333333333, 7.46, 7.453333333333333, 7.43351, 7.430000000000001, 7.419133333333334,
                           7.403333333333333, 7.4024366666666666, 7.400276666666667, 7.4, 7.4, 7.37, 7.36,
                           7.343333333333333, 7.33, 7.323333333333334, 7.322316666666667, 7.303333333333333,
                           7.300000000000001, 7.2981766666666665, 7.2968, 7.293333333333334, 7.275, 7.273333333333333,
                           7.2700000000000005, 7.266666666666667, 7.258990000000001, 7.256666666666667, 7.23, 7.23,
                           7.2299999999999995, 7.22, 7.22, 7.2122866666666665, 7.209999999999999, 7.203776666666666, 7.2,
                           7.2, 7.2, 7.196666666666667, 7.190083333333333, 7.175071666666666, 7.17, 7.169783333333334,
                           7.153333333333334, 7.12545, 7.122843333333334, 7.104733333333333, 7.1, 7.081026666666666, 7.07,
                           7.05198, 7.035706666666667, 7.03115, 7.03, 7.03, 7.03, 7.03, 7.03, 7.02, 7.013066666666667,
                           7.00767, 6.999553333333334, 6.99, 6.983333333333333, 6.963643333333334, 6.960000000000001,
                           6.95873, 6.951796666666667, 6.948796666666667, 6.948333333333333, 6.932113333333334,
                           6.913176666666667, 6.88, 6.88, 6.8527233333333335, 6.83, 6.83, 6.83, 6.826903333333333,
                           6.826666666666667, 6.808093333333334, 6.8, 6.791086666666667, 6.760000000000001, 6.75,
                           6.736666666666666, 6.730041666666667, 6.7299999999999995, 6.71898, 6.7113, 6.7, 6.7,
                           6.6913599999999995, 6.69, 6.681473333333334, 6.676806666666666, 6.658950000000001,
                           6.6543833333333335, 6.646666666666667, 6.640000000000001, 6.64, 6.62, 6.6160000000000005,
                           6.606865, 6.5933325, 6.58251, 6.5782799999999995, 6.57, 6.5600000000000005, 6.533966666666666,
                           6.53, 6.5266666666666655, 6.515406666666666, 6.510026666666667, 6.5, 6.489999999999999, 6.48406,
                           6.48, 6.48, 6.4799999999999995, 6.466666666666667, 6.458773333333334, 6.458333333333334, 6.44,
                           6.4, 6.4, 6.4, 6.4, 6.3999999999999995, 6.396626666666666, 6.383553333333333, 6.3781799999999995,
                           6.37, 6.36, 6.359999999999999, 6.359999999999999, 6.35873, 6.35479, 6.343333333333333,
                           6.333333333333333, 6.33244, 6.33, 6.33, 6.33, 6.31331, 6.311286666666667, 6.28,
                           6.276560000000001, 6.25, 6.239999999999999, 6.23, 6.226666666666667, 6.221733333333333, 6.22,
                           6.213333333333334, 6.21, 6.208333333333334, 6.200000000000001, 6.2, 6.197113333333333, 6.19344,
                           6.19, 6.178363333333333, 6.17034, 6.16, 6.1534466666666665, 6.145, 6.136746666666667,
                           6.129443333333334, 6.12, 6.093333333333334, 6.085193333333334, 6.0806, 6.073988333333334,
                           6.073333333333333, 6.070786666666667, 6.07, 6.0672066666666655, 6.065, 6.0600000000000005,
                           6.055413333333334, 6.05, 6.05, 6.033458333333333, 6.033333333333334, 6.03, 6.03, 6.03,
                           6.008579999999999, 6.008333333333334, 6.002776666666667, 6.0, 6.0, 5.97551, 5.96, 5.95,
                           5.9379800000000005, 5.930000000000001, 5.93, 5.926666666666667, 5.9, 5.9, 5.890000000000001,
                           5.883333333333334, 5.883006666666667, 5.859999999999999, 5.859999999999999, 5.85,
                           5.846666666666668, 5.845, 5.839203333333333, 5.837731666666667, 5.833333333333334,
                           5.833333333333334, 5.83, 5.811723333333333, 5.806900000000001, 5.805, 5.803333333333334, 5.8,
                           5.8, 5.785093333333333, 5.77, 5.760085, 5.759228333333334, 5.75589, 5.75, 5.74,
                           5.738476666666666, 5.73, 5.724003333333333, 5.723333333333334, 5.720400000000001,
                           5.706666666666667, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.7, 5.69, 5.69, 5.687508333333334, 5.677,
                           5.6725, 5.667050000000001, 5.663676666666667, 5.663333333333334, 5.6632766666666665,
                           5.6466666666666665, 5.6437, 5.643333333333334, 5.64, 5.634666666666666, 5.633333333333333,
                           5.631311666666667, 5.62, 5.609999999999999, 5.608996666666667, 5.6, 5.593333333333334, 5.59,
                           5.586666666666668, 5.576666666666667, 5.576666666666666, 5.57, 5.56, 5.54, 5.54, 5.53293,
                           5.531666666666667, 5.53, 5.526666666666667, 5.522606666666666, 5.520356666666666, 5.51982, 5.51,
                           5.508613333333334, 5.505703333333334, 5.5, 5.5, 5.5, 5.5, 5.473333333333334, 5.463519999999999,
                           5.456666666666667, 5.453333333333333, 5.449999999999999, 5.44, 5.433333333333334,
                           5.431996666666667, 5.425000000000001, 5.425000000000001, 5.415, 5.41, 5.403356666666667,
                           5.403333333333334, 5.4, 5.4, 5.4, 5.4, 5.4, 5.399333333333333, 5.390000000000001,
                           5.390000000000001, 5.380000000000001, 5.38, 5.38, 5.370733333333334, 5.3662399999999995,
                           5.351793333333333, 5.34, 5.33, 5.321253333333333, 5.32, 5.303333333333334, 5.303333333333334,
                           5.30116, 5.3, 5.296666666666667, 5.28482, 5.279999999999999, 5.270456666666666,
                           5.268460000000001, 5.265599999999999, 5.254886666666667, 5.25, 5.233333333333333,
                           5.216666666666666, 5.210593333333334, 5.2, 5.19381, 5.186666666666666, 5.183333333333334,
                           5.183333333333334, 5.18, 5.1771899999999995, 5.173333333333334, 5.17, 5.153333333333333, 5.15019,
                           5.1499999999999995, 5.149679999999999, 5.1486350000000005, 5.1450000000000005, 5.140000000000001,
                           5.14, 5.136666666666667, 5.130206666666667, 5.13, 5.1258333333333335, 5.124246666666667, 5.12,
                           5.118953333333334, 5.116666666666667, 5.109999999999999, 5.106666666666667, 5.085215000000001,
                           5.083333333333333, 5.076666666666667, 5.075920000000001, 5.056666666666667, 5.049376666666666,
                           5.0480599999999995, 5.045000000000001, 5.04, 5.03636, 5.035643333333334, 5.03, 5.03, 5.03, 5.03,
                           5.0240599999999995, 5.0200000000000005, 5.01, 5.003473333333334, 5.003333333333333, 5.0, 5.0,
                           5.0, 5.0, 5.0, 5.0, 5.0, 4.996666666666667, 4.99, 4.978333333333333, 4.963333333333334, 4.955,
                           4.950426666666667, 4.95, 4.95, 4.946666666666667, 4.944026666666667, 4.9383333333333335,
                           4.933333333333334, 4.928583333333334, 4.928, 4.923333333333334, 4.920000000000001, 4.92, 4.9,
                           4.9, 4.9, 4.9, 4.9, 4.893096666666667, 4.883333333333334, 4.876, 4.87, 4.8698500000000005,
                           4.866666666666667, 4.86059, 4.855026666666667, 4.8533333333333335, 4.84, 4.83872,
                           4.836666666666667, 4.835, 4.83, 4.8159, 4.813006666666666, 4.807556666666667, 4.800000000000001,
                           4.8, 4.8, 4.797000000000001, 4.795, 4.791989999999999, 4.790463333333333, 4.79, 4.78094,
                           4.777736666666668, 4.773486666666667, 4.773333333333333, 4.769048333333333, 4.765,
                           4.756666666666667, 4.75, 4.7448, 4.736773333333334, 4.73354, 4.73, 4.72927, 4.723333333333334,
                           4.716666666666667, 4.716666666666667, 4.7166266666666665, 4.713333333333334, 4.7114899999999995,
                           4.71, 4.706333333333333, 4.7, 4.7, 4.683333333333334, 4.683333333333334, 4.683303333333333,
                           4.674836666666667, 4.668685, 4.666666666666666, 4.663333333333333, 4.658333333333333,
                           4.6499999999999995, 4.649353333333334, 4.64247, 4.618953333333334, 4.611433333333333,
                           4.604627499999999, 4.5993, 4.596666666666667, 4.590380000000001, 4.573333333333334, 4.57144,
                           4.5600000000000005, 4.556666666666667, 4.5349450000000004, 4.534613333333334, 4.530766666666667,
                           4.529763333333333, 4.5207, 4.52, 4.511613333333333, 4.5, 4.479276666666667, 4.474666666666667,
                           4.461083333333334, 4.46, 4.456666666666667, 4.45, 4.435326666666667, 4.43, 4.423333333333334,
                           4.42, 4.416666666666667, 4.410888333333333, 4.40227, 4.4, 4.4, 4.3950000000000005,
                           4.390000000000001, 4.363333333333333, 4.3467400000000005, 4.338993333333333, 4.33, 4.33,
                           4.304556666666667, 4.301666666666667, 4.3, 4.3, 4.2825, 4.279613333333334, 4.266666666666667,
                           4.263333333333334, 4.261223333333333, 4.23064, 4.23, 4.23, 4.22816, 4.216666666666667,
                           4.203333333333333, 4.203333333333333, 4.201196666666667, 4.2, 4.2, 4.2, 4.1866666666666665,
                           4.167156666666666, 4.156666666666666, 4.140000000000001, 4.139943333333333, 4.136621666666667,
                           4.136183333333333, 4.13061, 4.12, 4.1136425, 4.1033333333333335, 4.098333333333333, 4.08, 4.08,
                           4.073333333333333, 4.0600000000000005, 4.043333333333333, 4.036666666666666, 4.023333333333333,
                           4.0, 4.0, 4.0, 3.9966666666666666, 3.9870000000000005, 3.98475, 3.9833333333333334,
                           3.9833333333333334, 3.9733333333333336, 3.9633333333333334, 3.96, 3.96, 3.945,
                           3.9424733333333335, 3.9333333333333336, 3.9269299999999996, 3.9215400000000002,
                           3.9084933333333334, 3.9, 3.9, 3.883186666666667, 3.8805333333333336, 3.856666666666667, 3.855,
                           3.8499999999999996, 3.836666666666667, 3.829276666666667, 3.82777, 3.816666666666667,
                           3.7918399999999997, 3.7913699999999997, 3.7842666666666664, 3.7750333333333335, 3.77, 3.74,
                           3.7333333333333334, 3.7, 3.6966666666666668, 3.6933333333333334, 3.6834033333333336,
                           3.6807999999999996, 3.667226666666667, 3.663333333333333, 3.65, 3.626, 3.624546666666667,
                           3.6199999999999997, 3.6, 3.5983866666666664, 3.5966666666666667, 3.5966666666666667,
                           3.5896966666666663, 3.589025, 3.5779750000000003, 3.5666516666666666, 3.5556433333333333, 3.54,
                           3.533333333333333, 3.5319866666666666, 3.5300000000000002, 3.5233333333333334,
                           3.5199999999999996, 3.5, 3.5, 3.4973033333333334, 3.493333333333333, 3.49, 3.486666666666667,
                           3.473333333333333, 3.4673516666666666, 3.4308166666666664, 3.4200000000000004, 3.42, 3.41587,
                           3.4033333333333333, 3.4, 3.3866666666666667, 3.3798866666666667, 3.3633333333333333, 3.33,
                           3.3066666666666666, 3.3033333333333332, 3.276666666666667, 3.2733333333333334,
                           3.2283333333333335, 3.162833333333334, 3.1575966666666666, 3.149716666666667, 3.131906666666667,
                           3.1249233333333333, 3.1118333333333332, 3.0966666666666662, 3.093333333333333, 3.083333333333333,
                           3.0627266666666664, 3.0595833333333333, 3.053053333333333, 3.043333333333333, 3.0275000000000003,
                           3.0033333333333334, 3.0, 2.9845333333333333, 2.977433333333333, 2.9699999999999998, 2.93, 2.9,
                           2.85, 2.85, 2.816666666666667, 2.810286666666667, 2.763333333333333, 2.7333333333333334,
                           2.7175000000000002, 2.71, 2.6975000000000002, 2.6900116666666665, 2.6280616666666665,
                           2.6090299999999997, 2.5978616666666667, 2.57, 2.5448983333333333, 2.54067, 2.52, 2.47795,
                           2.4443366666666666, 2.41, 2.3850000000000002, 2.27, 2.2441825, 2.2416666666666663,
                           2.1951799999999997, 2.1799999999999997, 2.163333333333333, 2.13, 2.1, 1.7993899999999998,
                           1.7650000000000001, 1.31, 4.3100233333333335, 10.2, 6.29243, 10.419293333333334,
                           5.168221666666668, 6.7, 7.74, 9.7, 8.03842, 7.73, 9.7, 5.674025, 6.06035, 4.298333333333333, 8.0,
                           3.589123333333333, 6.1439, 7.705913333333333, 11.39, 8.0, 3.8994133333333334, 4.37748,
                           4.353393333333334, 2.5772049999999997, 3.9458100000000003, 5.828693333333334, 3.9165, 2.05811,
                           4.929376666666666, 4.859116666666667, 5.5600000000000005, 8.03428, 6.047476666666666, 6.83,
                           8.23805, 6.46246, 8.541443333333333, 7.5279]
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