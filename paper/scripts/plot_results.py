import random

import pandas as pd
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
from scipy import stats
from statannot import add_stat_annotation
from AutoCaSc_core.tools import add_categories
import matplotlib.patches as mpatches
import numpy as np

# ROOT_DIR = "/home/johann/PycharmProjects/AutoCaSc_project_folder/update_data/data/"
# ROOT_DIR = "/Users/johannkaspar/OneDrive/Promotion/AutoCaSc_project_folder/update_data/data/"
ROOT_DIR = "/update_data/data/"

colors = ["#a8d0db", "#4062BB", "#FED766", "#a37a74"]
pal = sns.color_palette(colors)
sns.set_palette(pal)

sysid_primary_ensemble = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv", usecols=["Ensembl id"])["Ensembl id"].to_list()
sysid_candidates_ensemble = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv", usecols=["Ensembl id"])["Ensembl id"].to_list()
princeton_negative_ensemble = pd.read_csv(ROOT_DIR + "ASD_translated_to_ensembl.csv")["ensemble_id"].to_list()

sysid_primary_entrez = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv", usecols=["Entrez id"])["Entrez id"].to_list()
sysid_candidates_entrez = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv", usecols=["Entrez id"])["Entrez id"].to_list()
princeton_negative_entrez = pd.read_csv(ROOT_DIR + "ASD_translated_to_ensembl.csv")["entrez_id"].to_list()



global order
order = ["unknown", "negative control", "candidate", "known NDD"]

def select_validation_df(all_gene_data):
    sysid_primary = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv",
                                usecols=["Entrez id", "Ensembl id"])
    sysid_primary.columns = ["entrez_id", "ensemble_id"]
    sysid_candidates = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv",
                                   usecols=["Entrez id", "Ensembl id"])
    sysid_candidates.columns = ["entrez_id", "ensemble_id"]
    princeton_negative = pd.read_csv(ROOT_DIR + "ASD_translated_to_ensembl.csv")["entrez_id"].to_list()

    morbid_gene_symbols_list = pd.read_csv(ROOT_DIR + "MorbidGenes-Panel-v5_2020-08-26_for_varvis.csv",
                                           header=None).iloc[:, 0].to_list()
    all_genes_df = pd.read_csv(ROOT_DIR + "hgnc_protein_coding.tsv", index_col=False,
                               usecols=["entrez_id", "gene_symbol"], sep="\t",
                               dtype={"entrez_id": "Int32", "gene_symbol": str})
    morbid_genes = all_genes_df.loc[all_genes_df.gene_symbol.isin(morbid_gene_symbols_list)][
        "entrez_id"].dropna().to_list()
    panel_app_genes = pd.read_csv(ROOT_DIR + "Intellectual disability.tsv", sep="\t", usecols=["HGNC"])[
        "HGNC"].to_list()
    panel_app_genes = [int(x.strip("HGNC:")) for x in panel_app_genes if type(x) == str]
    g2p_dd_list = pd.read_csv(ROOT_DIR + "DDG2P_26_10_2020.csv",
                              usecols=["hgnc id"],
                              dtype={"hgnc id": "Int32"})
    g2p_dd_list = g2p_dd_list["hgnc id"].to_list()
    negative_gene_list = list(set(morbid_genes) - set(panel_app_genes)
                              - set(sysid_primary.entrez_id.to_list() + sysid_candidates.entrez_id.to_list())
                              - set(g2p_dd_list))

    random.seed(42)
    _negative_gene_list = random.sample(negative_gene_list, round(0.8 * len(negative_gene_list)))
    negative_control = list(set(negative_gene_list) - set(_negative_gene_list))
    _rows_id = random.sample(range(0, len(sysid_primary)), round(0.8 * len(sysid_primary)))
    rows_id = list(set(range(0, len(sysid_primary))) - set(_rows_id))
    sysid_primary = sysid_primary.loc[rows_id, :].reset_index(drop=True)
    _rows_id = random.sample(range(len(sysid_candidates)), round(0.8 * len(sysid_candidates)))
    rows_id = list(set(range(0, len(sysid_candidates))) - set(_rows_id))
    sysid_candidates = sysid_candidates.loc[rows_id, :].reset_index(drop=True)

    all_gene_data = all_gene_data.loc[all_gene_data.entrez_id.isin(sysid_primary.entrez_id.to_list()
                                                                   + sysid_candidates.entrez_id.to_list()
                                                                   + negative_control)]
    return all_gene_data

def plot_parameter_genescores(validation_run=False):
    global all_gene_data, order
    sns.set(style="ticks",
            # font_scale=0.4
            )
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(12, 8), sharey=True, sharex=True, dpi=300)
    all_gene_data = pd.read_csv(ROOT_DIR + "all_gene_data.csv")

    if validation_run:
        all_gene_data = select_validation_df(all_gene_data)
        order = ["negative control", "known NDD"]

    all_gene_data = add_categories(all_gene_data, "entrez_id", "entrez", "morbid_genes")
    # gtex_data = gtex_data.loc[:4000]

    # axs[0, 0] = get_parameter_plot("gtex_score")
    # axs[0, 1] = get_parameter_plot("disgenet_score")
    # axs[0, 2] = get_parameter_plot("denovo_rank_score")
    # axs[1, 0] = get_parameter_plot("mgi_score")
    # axs[1, 1] = get_parameter_plot("pubtator_score")
    # axs[1, 2] = get_parameter_plot("string_score")
    n = 0
    parameters = {"gtex_score": "GTEx score",
                  "disgenet_score": "DisGeNET score",
                  "denovo_rank_score": "PsyMuKB rank score",
                  "mgi_score": "MGI score",
                  "pubtator_score": "PubTator Central score",
                  "string_score": "StringDB score"}

    for n, _parameter in enumerate(parameters.keys()):
        row = n//3
        col = n%3
        subplot = sns.stripplot(x="sys_category", y=_parameter, data=all_gene_data,
                       palette=pal,
                       size=2,
                       alpha=0.4,
                       order=order,
                       zorder=-10,
                       jitter=0.35,
                       ax=axs[row][col]
                           )
        subplot = sns.boxplot(x="sys_category", y=_parameter, data=all_gene_data,
                     showcaps=True,
                     showbox=False,
                     boxprops={'facecolor': '#D3D3D3',
                               "color":"#D3D3D3"},
                     whis=0,
                     showfliers=False,
                     order=order,
                     width=0.1,
                     zorder=5,
                      ax=axs[row][col]
                     )
        subplot = sns.boxplot(x="sys_category", y=_parameter, data=all_gene_data,
                         showcaps=False,
                         showbox=True,
                         # boxprops={'facecolor': 'None'},
                         whis=0,
                         showfliers=False,
                         whiskerprops={'linewidth': 0.01},
                         order=order,
                         width=0.01,
                         zorder=2.6,
                      ax=axs[row][col]
                         )
        subplot.set_xticklabels(subplot.get_xticklabels(), rotation=45, horizontalalignment='right')
        subplot.set(xlabel=None, ylabel=None, title=parameters.get(_parameter))

    fig.suptitle('Genescores by gene category', fontsize=10)
    plt.show()

def weighted_score(validation_run=False):
    all_gene_data = pd.read_csv(ROOT_DIR + "all_gene_data.csv")
    all_gene_data = all_gene_data.drop_duplicates(subset=["entrez_id"])

    if validation_run:
        all_gene_data = select_validation_df(all_gene_data)
        order = ["negative control", "known NDD"]

    all_gene_data = add_categories(all_gene_data, "entrez_id", "entrez", "morbid_genes")

    # ax = sns.violinplot(x="sys_category", y="weighted_score", data=all_gene_data,
    #                  order=order, color="0.5")
    # ax = sns.stripplot(x="sys_category", y="weighted_score", data=all_gene_data, jitter=True,
    #                  order=order, alpha=0.2)

    ax = plt.figure(figsize=(4.5, 3))

    # scatter
    ax = sns.stripplot(x="sys_category", y="weighted_score", data=all_gene_data,
                            palette=pal,
                            size=5,
                            alpha=0.6,
                            order=order,
                            zorder=-10,
                            jitter=0.2
                            )
    # whiskers
    ax = sns.boxplot(x="sys_category", y="weighted_score", data=all_gene_data,
                          showcaps=True,
                          showbox=False,
                          boxprops={'facecolor': '#D3D3D3',
                                    "color": "#D3D3D3"},
                          whis=0,
                          showfliers=False,
                          order=order,
                          width=0.2,
                          zorder=5
                          )

    # boxline
    ax = sns.boxplot(x="sys_category", y="weighted_score", data=all_gene_data,
                          showcaps=False,
                          showbox=True,
                          whis=0,
                          showfliers=False,
                          whiskerprops={'linewidth': 0},
                          order=order,
                          width=0.01,
                          zorder=2.6
                          )

    # significance
    add_stat_annotation(ax, data=all_gene_data, x="sys_category", y="weighted_score", order=order,
                        box_pairs=[("negative control", "known NDD")],
                        test='Mann-Whitney', text_format='star', loc='outside', line_offset_to_box=0.001,
                        line_height=0.05, text_offset=2, verbose=2)

    ax.set(ylim=(0, 1.4))

    ax.set_title("Genescores by SysID category")
    ax.set_ylabel("Weighted score")
    ax.yaxis.set_ticks(np.arange(0.0, 1.1, 0.2))
    ax.set_xlabel("")
    plt.xticks(color='w')
    ax.xaxis.set_ticks([])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # set legend
    negative_patch = mpatches.Patch(color=colors[0], label='Negative control')
    positive_patch = mpatches.Patch(color=colors[1], label='NDD control')
    plt.legend(handles=[negative_patch, positive_patch],
               framealpha=0.5,
               fontsize="small",
               bbox_to_anchor=(1.05, 1.0))
    ax.get_figure().savefig("/Users/johannkaspar/OneDrive/Promotion/AutoCaSc_project_folder"
                            "/Manuskript/items/Fig 2/weighted_sum_plot_seed_42.png")
    plt.show()

def plot_candidate_scores():
    candidate_scores = pd.read_csv("/Users/johannkaspar/OneDrive/Promotion/AutoCaSc_project_folder/paper/data/manual_vs_automated.csv",
                                   sep=";")
    candidate_scores = candidate_scores.loc[candidate_scores.AutoCaSc != "-"]
    candidate_scores.loc[:, "CaSc"] = candidate_scores.loc[:, "CaSc"].apply(lambda x: float(str(x).replace(",", ".")))
    candidate_scores.loc[:, "AutoCaSc"] = candidate_scores.loc[:, "AutoCaSc"].apply(lambda x: float(str(x).replace(",", ".")))

    plt.rcParams.update({'figure.figsize': (6, 6)})

    sns.lmplot(x='CaSc', y='AutoCaSc', data=candidate_scores,
               scatter_kws={'alpha':0.3},
               line_kws={"color":"black",
                         "alpha":0.7},
               fit_reg=True)
    plt.title("Correlation CaSc vs. AutoCaSc")
    plt.xlabel("CaSc")
    plt.ylabel("AutoCaSc")

    plt.xlim(0,16)
    plt.ylim(0,16)
    plt.xticks(range(16), [str(x) for x in range(16)])
    plt.yticks(range(16), [str(x) for x in range(16)])

    # plt.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    r, p = scipy.stats.pearsonr(candidate_scores.CaSc, candidate_scores.AutoCaSc)
    print(f"r = {r}\tp = {p}")
    r_patch = mpatches.Patch(color="white", label=f"r^2 = {round(r**2, 3)}")
    base = str(round(float(str(p).split("e")[0]), 3))
    exponent = str(p).split("e")[1]
    p_short = f"{base}e{exponent}"
    p_patch = mpatches.Patch(color="white", label=f"p = {p_short}")
    plt.legend(handles=[r_patch, p_patch],
               handlelength=0.0)
    plt.savefig("/Users/johannkaspar/OneDrive/Promotion/AutoCaSc_project_folder/Manuskript/Fig 2/correlation_plot.svg")
    plt.show()

def plot_pubtator_clean():
    sysid_primary = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv", usecols=["Entrez id"])["Entrez id"].to_list()
    sysid_candidates = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv", usecols=["Entrez id"])[
        "Entrez id"].to_list()
    princeton_negative = pd.read_csv(ROOT_DIR + "ASD_translated_to_ensembl.csv")["gene id"].to_list()

    pubtator = pd.read_csv(ROOT_DIR + "pubtator_central/gene_scores/gene_scores_p_cutoff_0,0001_clean.csv")

    pubtator["sys_primary"] = pubtator.gene_id.isin(sysid_primary).astype(int)
    pubtator["sys_candidate"] = pubtator.gene_id.isin(sysid_candidates).astype(int)
    pubtator["sys"] = pubtator.gene_id.isin(sysid_primary + sysid_candidates).astype(int)
    pubtator["sys_category"] = "unknown"
    pubtator.loc[pubtator.sys_candidate == 1, "sys_category"] = "candidate"
    pubtator.loc[pubtator.sys_primary == 1, "sys_category"] = "known NDD"
    pubtator.loc[pubtator.gene_id.isin(princeton_negative), "sys_category"] = "negative control"

    order = ["unknown", "negative control", "candidate", "known NDD"]

    ax = plt.figure(figsize=(6,6))
    ax = sns.boxplot(x="sys_category", y="gene_score", data=pubtator, showfliers=False,
                     order=order)
    add_stat_annotation(ax, data=pubtator, x="sys_category", y="gene_score", order=order,
                       box_pairs=[("unknown", "candidate"), ("candidate", "negative control"), ("candidate", "known NDD")],
                       test='Mann-Whitney', text_format='simple', loc='outside', line_offset_to_box=0.001,
                        line_height=0.05, text_offset=2, verbose=2)
    ax.set(ylim=(0, 1800))
    # ax.set(ylim=(0, 0.07))

    ax.set_title(f"pubtator gene scores")
    ax.set_xlabel("SysID category")
    ax.set_ylabel(f"gene score")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_figure().savefig(ROOT_DIR + f"pubtator_central/plot_pubtator_clean.png")
    plt.show()

# plot_gtex_brain_sum()
# plot_string()
# plot_psymukb()
# plot_pubtator()
# plot_string()
# plot_disgenet()
# plot_mgi()
# rain_disgenet()
# weighted_score(validation_run=True)
# plot_candidate_scores()

plot_parameter_genescores(validation_run=False)