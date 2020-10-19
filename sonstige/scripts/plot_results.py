import pandas as pd
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
from scipy import stats
from statannot import add_stat_annotation
from tools import add_categories
import ptitprince as pt
import matplotlib.patches as mpatches
from sklearn.metrics import r2_score

pal = sns.color_palette(["#a8d0db", "#4062BB", "#FED766", "#a37a74"])
sns.set_palette(pal)


# old functions
def plot_string_sys_sum(parameter="sys_sum"):
    string_processed = pd.read_csv("/home/johann/AutoCaSc/data/string/interaction_sum_sys_all.csv")
    string_processed["sys"] = string_processed.gene_id.isin(sysid_candidates + sysid_primary).astype(int)
    string_processed["sys_primary"] = string_processed.gene_id.isin(sysid_primary).astype(int)
    string_processed["sys_candidate"] = string_processed.gene_id.isin(sysid_candidates).astype(int)
    string_processed["negative"] = string_processed.gene_id.isin(princeton_negative).astype(int)
    string_processed["sys_category"] = "unknown"
    string_processed.loc[string_processed.sys_candidate == 1, "sys_category"] = "candidate"
    string_processed.loc[string_processed.sys_primary == 1, "sys_category"] = "known NDD"
    string_processed.loc[string_processed.negative == 1, "sys_category"] = "negative control"


    # for parameter in list(string_processed.columns)[1:7]:
    #     print(stats.mannwhitneyu(string_processed.loc[string_processed.sys_candidate == 1][parameter].to_list(),
    #     string_processed.loc[string_processed.sys == 0][parameter].to_list(), alternative="greater"))

    order = ["unknown", "negative control", "candidate", "known NDD"]

    ax = plt.figure(figsize=(7,9))
    ax = sns.boxplot(x="sys_category", y=parameter, data=string_processed, showfliers=False,
                     order=order)
    add_stat_annotation(ax, data=string_processed, x="sys_category", y=parameter, order=order,
                       box_pairs=[("unknown", "candidate"), ("unknown", "known NDD"), ("candidate", "known NDD")],
                       test='Mann-Whitney', text_format='full', loc='outside', line_offset_to_box=0.001,
                        line_height=0.05, text_offset=2, verbose=2)
    ax.set(ylim=(0, 150))

    ax.set_title(f"sum of interactions with SysID genes in StringDB")
    ax.set_xlabel("SysID category")
    ax.set_ylabel(f"corrected sum of interaction confidences with SysID genes")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_figure().savefig(ROOT_DIR + f"string/plot_{parameter}.png")
    plt.show()

def plot_gtex_sum_hue():
    colors = ["#fadf63", "#3877B6"]
    gtex_data = pd.read_csv("/home/johann/AutoCaSc/data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                                skiprows=2)

    gtex_data["ensemble_id"] = gtex_data["Name"].str.split(".")[0]
    gtex_data = gtex_data.drop(columns="Name")

    # gtex_data = gtex_data.loc[gtex_data.brain_sum != 0]
    gtex_data["sys_category"] = "unknown"
    gtex_data.loc[gtex_data.sys_candidate == 1, "sys_category"] = "candidate"
    gtex_data.loc[gtex_data.sys_primary == 1, "sys_category"] = "known NDD"

    gtex_temp_brain = gtex_data[["ensemble_id", "sys_primary", "sys_candidate", "sys", "brain_sum", "ratio", "brain_max", "brain_median", "sys_category"]]
    gtex_temp_brain["tissue"] = "brain"
    gtex_temp_brain.columns = ["ensemble_id", "sys_primary", "sys_candidate", "sys", "sum", "ratio", "max", "median", "sys_category", "tissue"]

    gtex_temp_other = gtex_data[["ensemble_id", "sys_primary", "sys_candidate", "sys", "other_sum", "ratio", "other_max", "other_median", "sys_category"]]
    gtex_temp_other["tissue"] = "other"
    gtex_temp_other.columns = ["ensemble_id", "sys_primary", "sys_candidate", "sys", "sum", "ratio", "max", "median", "sys_category", "tissue"]

    gtex_temp = pd.concat([gtex_temp_brain, gtex_temp_other], ignore_index=True)

    sns.set_palette(sns.color_palette(colors))
    ax = sns.boxplot(x="sys_category", y="sum", hue="tissue", data=gtex_temp, showfliers=False, order=["unknown", "candidate", "known NDD"])

    ax.set_title("tissue expression by SysID category")
    ax.set_xlabel("SysID category")
    ax.set_ylabel("Sum of median TPM")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_figure().savefig(ROOT_DIR + "gtex/plot_sum_hue.png")
    plt.show()

def plot_gtex_ratio(parameter="ratio"):
    colors = ["windows blue", "amber", "faded green", "greyish", "dusty purple"]
    sns.set_palette(sns.xkcd_palette(colors))

    gtex_data = pd.read_csv("/home/johann/AutoCaSc/data/gtex/gtex_data.csv")

    # gtex_data = gtex_data.loc[gtex_data.brain_sum != 0]
    gtex_data["sys_category"] = "unknown"
    gtex_data.loc[gtex_data.sys_candidate == 1, "sys_category"] = "candidate"
    gtex_data.loc[gtex_data.sys_primary == 1, "sys_category"] = "known NDD"
    gtex_data.loc[gtex_data.ensemble_id.isin(princeton_negative_ensemble), "sys_category"] = "negative control"


    ax = plt.figure(figsize=(6,6))
    ax = sns.boxplot(x="sys_category", y=parameter, data=gtex_data, showfliers=False,
                     order=order)

    add_stat_annotation(ax, data=gtex_data, x="sys_category", y=parameter, order=order,
                        box_pairs=[("unknown", "candidate"), ("unknown", "known NDD"),
                                   ("negative control", "candidate")],
                        test='Mann-Whitney', text_format='simple', loc='outside', line_offset_to_box=0.001,
                        line_height=0.05, text_offset=2, verbose=2)
    ax.set(ylim=(0, 1.2))

    ax.set_title(f"{parameter} expression by SysID category")
    ax.set_xlabel("SysID category")
    ax.set_ylabel(f"{parameter} of median TPM")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_figure().savefig(ROOT_DIR + f"gtex/plot_{parameter}.png")
    plt.show()

def plot_gtex_all_sum(parameter="all_sum"):
    colors = ["windows blue", "amber", "faded green", "greyish", "dusty purple"]
    sns.set_palette(sns.xkcd_palette(colors))

    gtex_data = pd.read_csv("/home/johann/AutoCaSc/data/gtex/gtex_data.csv")

    # gtex_data = gtex_data.loc[gtex_data.brain_sum != 0]
    gtex_data["sys_category"] = "unknown"
    gtex_data.loc[gtex_data.sys_candidate == 1, "sys_category"] = "candidate"
    gtex_data.loc[gtex_data.sys_primary == 1, "sys_category"] = "known NDD"
    gtex_data.loc[gtex_data.ensemble_id.isin(princeton_negative_ensemble), "sys_category"] = "negative control"


    ax = plt.figure(figsize=(6,7))
    ax = sns.boxplot(x="sys_category", y=parameter, data=gtex_data, showfliers=False,
                     order=order)

    add_stat_annotation(ax, data=gtex_data, x="sys_category", y=parameter, order=order,
                        box_pairs=[("unknown", "candidate"), ("unknown", "known NDD"), ("candidate", "known NDD")],
                        test='Mann-Whitney', text_format='full', loc='outside', line_offset_to_box=0.001,
                        line_height=0.05, text_offset=2, verbose=2)
    ax.set(ylim=(0, 5500))

    ax.set_title(f"{parameter} expression by SysID category")
    ax.set_xlabel("SysID category")
    ax.set_ylabel(f"{parameter} of median TPM")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_figure().savefig(ROOT_DIR + f"gtex/plot_{parameter}.png")
    plt.show()


# global ROOT_DIR
order = ["unknown", "negative control", "candidate", "known NDD"]

ROOT_DIR = "/Users/johannkaspar/OneDrive/Promotion/AutoCaSc_project_folder/AutoCaSc_maintenance/data/"
# ROOT_DIR = "/home/johann/AutoCaSc/data/"

sysid_primary_ensemble = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv", usecols=["Ensembl id"])["Ensembl id"].to_list()
sysid_candidates_ensemble = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv", usecols=["Ensembl id"])["Ensembl id"].to_list()
princeton_negative_ensemble = pd.read_csv(ROOT_DIR + "ASD_translated_to_ensembl.csv")["Gene stable ID"].to_list()

sysid_primary_entrez = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv", usecols=["Entrez id"])["Entrez id"].to_list()
sysid_candidates_entrez = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv", usecols=["Entrez id"])["Entrez id"].to_list()
princeton_negative_entrez = pd.read_csv(ROOT_DIR + "ASD_translated_to_ensembl.csv")["gene id"].to_list()

def plot_psymukb():
    denovo_df = pd.read_csv("/home/johann/AutoCaSc/data/psymukb/all_genes_denovo_NDD.csv")
    denovo_df = add_categories(denovo_df, "entrez_id", "entrez")

    ax = plt.figure(figsize=(6, 6))
    # ax = sns.boxplot(x="sys_category", y="denovo_count", data=denovo_df, showfliers=False,
    #                  order=order)
    print("1")
    ax = sns.swarmplot(x="sys_category", y="denovo_count", data=denovo_df, order=order)
    print("2")
    # add_stat_annotation(ax, data=denovo_df, x="sys_category", y="denovo_count", order=order,
    #                     box_pairs=[("unknown", "candidate"), ("unknown", "known NDD"),
    #                                ("negative control", "candidate")],
    #                     test='Mann-Whitney', text_format='simple', loc='outside', line_offset_to_box=0.001,
    #                     line_height=0.05, text_offset=2, verbose=2)
    plt.show()


def plot_gtex_brain_sum(parameter="brain_sum"):
    # colors = ["windows blue", "amber", "faded green", "greyish", "dusty purple"]
    # sns.set_palette(sns.xkcd_palette(colors))
    sns.set_style("ticks")

    gtex_data = pd.read_csv(ROOT_DIR + "gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
        skiprows=2,
        sep="\t")
    gtex_data["ensemble_id"] = gtex_data["Name"].apply(lambda x: x.split(".")[0])
    gtex_data = gtex_data.drop(columns="Name")

    columns = gtex_data.columns
    brain_columns = [x for x in columns if "Brain" in x]
    gtex_data["brain_sum"] = gtex_data[brain_columns].sum(axis=1)

    gtex_data = gtex_data.loc[gtex_data.brain_sum != 0]
    gtex_data = add_categories(gtex_data, "ensemble_id", "ensemble")

    ax = plt.figure(figsize=(6,4))
    ax = sns.boxplot(x="sys_category", y=parameter, data=gtex_data, showfliers=False,
                     order=order)
    # sns.despine(offset=10, trim=True);


    add_stat_annotation(ax, data=gtex_data, x="sys_category", y=parameter, order=order,
                        box_pairs=[("unknown", "candidate"), ("unknown", "known NDD"),
                                   ("negative control", "candidate")],
                        test='Mann-Whitney', text_format='simple', loc='outside', line_offset_to_box=0.0001,
                        line_height=0.02, text_offset=1.2, verbose=2)
    ax.set(ylim=(0, 1100))

    ax.set_title(f"{parameter} expression by SysID category")
    ax.set_xlabel("SysID category")
    ax.set_ylabel(f"{parameter} of median TPM")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.get_figure().savefig("/home/johann/AutoCaSc/manuscript/Fig 2/plots/gtex.png")
    plt.show()

def plot_string(parameter="sum_corrected_exp_07"):
    string_processed = pd.read_csv("/home/johann/AutoCaSc/data/string/gene_scores.csv")
    string_processed = add_categories(string_processed, "gene_id", "ensemble")

    ax = plt.figure(figsize=(7,9))
    ax = sns.boxplot(x="sys_category", y=parameter, data=string_processed, showfliers=False,
                     order=order)
    add_stat_annotation(ax, data=string_processed, x="sys_category", y=parameter, order=order,
                       box_pairs=[("unknown", "candidate"), ("unknown", "known NDD"), ("negative control", "candidate")],
                       test='Mann-Whitney', text_format='full', loc='outside', line_offset_to_box=0.001,
                        line_height=0.05, text_offset=2, verbose=2)
    ax.set(ylim=(0, 1300))

    ax.set_xlabel("SysID category")
    ax.set_ylabel(f"corrected sum of interaction confidences with SysID genes")
    ax.set_title(f"corrected sum of interactions with SysID genes in StringDB")


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_figure().savefig("/home/johann/AutoCaSc/manuscript/Fig 2/plots/string.png")
    plt.show()

def plot_disgenet():
    disgenet = pd.read_csv("/home/johann/AutoCaSc/data/disgenet/disgenet_gene_scores.csv")
    disgenet = add_categories(disgenet, "entrez_id", "entrez")

    ax = plt.figure(figsize=(6, 6))
    ax = sns.boxplot(x="sys_category", y="gene_score", data=disgenet, showfliers=False,
                     order=order)
    add_stat_annotation(ax, data=disgenet, x="sys_category", y="gene_score", order=order,
                        box_pairs=[("unknown", "candidate"), ("unknown", "known NDD"), ("candidate", "known NDD")],
                        test='Mann-Whitney', text_format='full', loc='outside', line_offset_to_box=0.001,
                        line_height=0.05, text_offset=2, verbose=2)
    ax.set(ylim=(0, 1.5))

    ax.set_title(f"disgenet gene scores")
    ax.set_xlabel("SysID category")
    ax.set_ylabel(f"gene score")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_figure().savefig("/home/johann/AutoCaSc/manuscript/Fig 2/plots/disgenet.png")
    plt.show()

def rain_disgenet():
    sns.set(style="whitegrid", font_scale=1.5)
    disgenet = pd.read_csv(ROOT_DIR + "disgenet/disgenet_gene_scores.csv")
    disgenet = add_categories(disgenet, "entrez_id", "entrez")

    f, ax = plt.subplots(figsize=(8, 12))
    dy = "sys_category"; dx="gene_score"; ort="h"; pal="Set2"; sigma=0.1

    ax = pt.half_violinplot(x=dx, y=dy, data=disgenet, palette=pal,
                            bw=sigma, cut=0., scale="area", width=3, orient=ort, inner=None)
    ax = sns.stripplot(x=dx, y=dy, data=disgenet, palette=pal, edgecolor="white", size=5,
                       jitter=1, zorder=0, orient=ort, alpha=0.3)
    ax = sns.boxplot(x=dx, y=dy, data=disgenet, color="black", width=0.15, zorder=10, showcaps=True,
                     boxprops={"facecolor":"none", "zorder":10}, showfliers=False,
                     whiskerprops={"linewidth":2, "zorder":10}, saturation=1, orient=ort)

    # ax = plt.figure(figsize=(6, 6))
    # ax = sns.boxplot(x="sys_category", y="gene_score", data=disgenet, showfliers=False,
    #                  order=order)
    # add_stat_annotation(ax, data=disgenet, x="sys_category", y="gene_score", order=order,
    #                     box_pairs=[("unknown", "candidate"), ("unknown", "known NDD"), ("candidate", "known NDD")],
    #                     test='Mann-Whitney', text_format='full', loc='outside', line_offset_to_box=0.001,
    #                     line_height=0.05, text_offset=2, verbose=2)
    # ax.set(ylim=(0, 1.5))
    #
    # ax.set_title(f"disgenet gene scores")
    # ax.set_xlabel("SysID category")
    # ax.set_ylabel(f"gene score")
    #
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.get_figure().savefig("/home/johann/AutoCaSc/manuscript/Fig 2/plots/disgenet.png")
    plt.show()

def plot_parameter_genescores():
    global all_gene_data
    sns.set(style="ticks", font_scale=1)
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(10, 10), sharey=True, sharex=True, dpi=300)
    all_gene_data = pd.read_csv(ROOT_DIR + "all_gene_data.csv")
    all_gene_data = add_categories(all_gene_data, "entrez_id", "entrez")
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
                     width=0.2,
                     zorder=5,
                      ax=axs[row][col]
                     )
        subplot = sns.boxplot(x="sys_category", y=_parameter, data=all_gene_data,
                         showcaps=False,
                         showbox=True,
                         # boxprops={'facecolor': 'None'},
                         whis=0,
                         showfliers=False,
                         whiskerprops={'linewidth': 0},
                         order=order,
                         width=0.01,
                         zorder=2.6,
                      ax=axs[row][col]
                         )
        subplot.set_xticklabels(subplot.get_xticklabels(), rotation=45, horizontalalignment='right')
        subplot.set(xlabel=None, ylabel=None, title=parameters.get(_parameter))

    fig.suptitle('Genescores by gene category', fontsize=20)

    # ---------------
    # dy="sys_category"; dx="gtex_score"; ort="h"; sigma=0.05
    # ax = pt.half_violinplot(x=dx, y=dy, data=gtex_data, palette=pal,
    #                         bw=sigma, cut=0., scale="area", width=2.5, orient=ort, inner=None)
    # ax = sns.stripplot(x=dx, y=dy, data=gtex_data, palette=pal, edgecolor="white", size=5,
    #                    jitter=1, zorder=0, orient=ort, alpha=0.1)
    # ax = sns.boxplot(x=dx, y=dy, data=gtex_data, color="black", width=0.15, zorder=10, showcaps=True,
    #                  boxprops={"facecolor":"none", "zorder":10}, showfliers=False,
    #                  whiskerprops={"linewidth":2, "zorder":10}, saturation=1, orient=ort)

    # ax = plt.figure(figsize=(6, 6))
    # ax = sns.boxplot(x="sys_category", y="gene_score", data=disgenet, showfliers=False,
    #                  order=order)
    # add_stat_annotation(ax, data=disgenet, x="sys_category", y="gene_score", order=order,
    #                     box_pairs=[("unknown", "candidate"), ("unknown", "known NDD"), ("candidate", "known NDD")],
    #                     test='Mann-Whitney', text_format='full', loc='outside', line_offset_to_box=0.001,
    #                     line_height=0.05, text_offset=2, verbose=2)
    # ax.set(ylim=(0, 1.5))
    #
    # ax.set_title(f"disgenet gene scores")
    # ax.set_xlabel("SysID category")
    # ax.set_ylabel(f"gene score")
    #
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.get_figure().savefig("/home/johann/AutoCaSc/manuscript/Fig 2/plots/disgenet.png")
    plt.savefig("/Users/johannkaspar/OneDrive/Promotion/AutoCaSc_project_folder/Manuskript/Fig 2/gene_scores.png")
    plt.show()


def plot_mgi():
    mgi_df = pd.read_csv("/home/johann/AutoCaSc/data/mgi/mgi_gene_scores.csv")
    mgi_df = add_categories(mgi_df, "entrez_id", "entrez")

    ax = plt.figure(figsize=(6, 6))
    ax = sns.boxplot(x="sys_category", y="gene_score", data=mgi_df, showfliers=False,
                     order=order)
    add_stat_annotation(ax, data=mgi_df, x="sys_category", y="gene_score", order=order,
                        box_pairs=[("unknown", "candidate"), ("unknown", "known NDD"), ("candidate", "known NDD")],
                        test='Mann-Whitney', text_format='full', loc='outside', line_offset_to_box=0.001,
                        line_height=0.05, text_offset=2, verbose=2)
    ax.set(ylim=(0, 7))

    ax.set_title(f"MGI gene score")
    ax.set_xlabel("SysID category")
    ax.set_ylabel(f"gene score")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_figure().savefig("/home/johann/AutoCaSc/manuscript/Fig 2/plots/mgi.png")
    plt.show()



# def plot_pubtator():
#     sysid_primary = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv", usecols=["Entrez id"])["Entrez id"].to_list()
#     sysid_candidates = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv", usecols=["Entrez id"])[
#         "Entrez id"].to_list()
#     princeton_negative = pd.read_csv(ROOT_DIR + "ASD_translated_to_ensembl.csv")["gene id"].to_list()
#
#     pubtator = pd.read_csv("/home/johann/AutoCaSc/data/pubtator_central/gene_ranks.csv")
#
#     ###################
#     # pubtator = pd.read_csv("/home/johann/AutoCaSc/data/pubtator_central/gene_scores/gene_scores_p_cutoff_0,0001_exp_1_pmid_single.csv")
#     ###################
#
#
#     pubtator["sys_primary"] = pubtator.gene_id.isin(sysid_primary).astype(int)
#     pubtator["sys_candidate"] = pubtator.gene_id.isin(sysid_candidates).astype(int)
#     pubtator["sys"] = pubtator.gene_id.isin(sysid_primary + sysid_candidates).astype(int)
#     pubtator["sys_category"] = "unknown"
#     pubtator.loc[pubtator.sys_candidate == 1, "sys_category"] = "candidate"
#     pubtator.loc[pubtator.sys_primary == 1, "sys_category"] = "known NDD"
#     pubtator.loc[pubtator.gene_id.isin(princeton_negative), "sys_category"] = "negative control"
#
#     order = ["unknown", "negative control", "candidate", "known NDD"]
#
#     ax = plt.figure(figsize=(6,6))
#     ax = sns.boxplot(x="sys_category", y="gene_score", data=pubtator, showfliers=False,
#                      order=order)
#     add_stat_annotation(ax, data=pubtator, x="sys_category", y="gene_score", order=order,
#                        box_pairs=[("unknown", "candidate"), ("candidate", "negative control"), ("candidate", "known NDD")],
#                        test='Mann-Whitney', text_format='simple', loc='outside', line_offset_to_box=0.001,
#                         line_height=0.05, text_offset=2, verbose=2)
#     ax.set(ylim=(0, 12))
#     # ax.set(ylim=(0, 0.07))
#
#     ax.set_title(f"pubtator gene scores")
#     ax.set_xlabel("SysID category")
#     ax.set_ylabel(f"gene score")
#
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#     ax.get_figure().savefig(ROOT_DIR + f"pubtator_central/plot_pubtator.png")
#     plt.show()
def plot_pubtator():
    pubtator_df = pd.read_csv("/home/johann/AutoCaSc/data/pubtator_central/gene_scores.csv")
    pubtator_df = add_categories(pubtator_df, "entrez_id", "entrez")

    ax = plt.figure(figsize=(6,6))
    ax = sns.boxplot(x="sys_category", y="pubtator_score", data=pubtator_df, showfliers=False,
                     order=order)
    add_stat_annotation(ax, data=pubtator_df, x="sys_category", y="pubtator_score", order=order,
                       box_pairs=[("unknown", "candidate"), ("candidate", "negative control"), ("candidate", "known NDD")],
                       test='Mann-Whitney', text_format='simple', loc='outside', line_offset_to_box=0.001,
                        line_height=0.05, text_offset=2, verbose=2)
    ax.set(ylim=(0, 50))
    # ax.set(ylim=(0, 0.07))

    ax.set_title(f"pubtator gene scores")
    ax.set_xlabel("SysID category")
    ax.set_ylabel(f"gene score")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_figure().savefig("/home/johann/AutoCaSc/manuscript/Fig 2/plots/pubtator.png")
    plt.show()

def weighted_score():
    all_gene_data = pd.read_csv("/home/johann/AutoCaSc/data/all_gene_data.csv")
    all_gene_data = all_gene_data.drop_duplicates(subset=["entrez_id"])

    all_gene_data.loc[all_gene_data.sys_category == "negative_control", "sys_category"] = "unknown"

    # order = ["unknown", "negative control", "candidate", "known NDD"]
    order = ["unknown", "candidate", "known NDD"]

    # ax = plt.figure(figsize=(3,6))
    ax = sns.violinplot(x="sys_category", y="weighted_score", data=all_gene_data,
                     order=order, color="0.5")
    ax = sns.stripplot(x="sys_category", y="weighted_score", data=all_gene_data, jitter=True,
                     order=order, alpha=0.2)
    add_stat_annotation(ax, data=all_gene_data, x="sys_category", y="weighted_score", order=order,
                       box_pairs=[("unknown", "candidate"), ("candidate", "known NDD")],
                       test='Mann-Whitney', text_format='star', loc='outside', line_offset_to_box=0.001,
                        line_height=0.05, text_offset=2, verbose=2)
    ax.set(ylim=(0, 1.5))

    ax.set_title("Genescores by SysID category")
    ax.set_xlabel("SysID category")
    ax.set_ylabel("Weighted score")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_figure().savefig(ROOT_DIR + "weighted_sum_plot.png")
    plt.show()

def plot_candidate_scores():
    candidate_scores = pd.read_csv("/Users/johannkaspar/OneDrive/Promotion/AutoCaSc_project_folder/sonstige/data/manual_vs_automated.csv",
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
    # print(r2_score(candidate_scores.CaSc, candidate_scores.AutoCaSc))
    plt.savefig("/Users/johannkaspar/OneDrive/Promotion/AutoCaSc_project_folder/Manuskript/Fig 2/correlation_plot.svg")
    plt.show()

def plot_pubtator_clean():
    sysid_primary = pd.read_csv(ROOT_DIR + "sysid/sysid_primary.csv", usecols=["Entrez id"])["Entrez id"].to_list()
    sysid_candidates = pd.read_csv(ROOT_DIR + "sysid/sysid_candidates.csv", usecols=["Entrez id"])[
        "Entrez id"].to_list()
    princeton_negative = pd.read_csv(ROOT_DIR + "ASD_translated_to_ensembl.csv")["gene id"].to_list()

    pubtator = pd.read_csv("/home/johann/AutoCaSc/data/pubtator_central/gene_scores/gene_scores_p_cutoff_0,0001_clean.csv")

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
# weighted_score()
# plot_candidate_scores()
plot_parameter_genescores()