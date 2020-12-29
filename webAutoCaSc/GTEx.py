import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import requests
from AutoCaSc_core.tools import safe_get
import time

import pandas as pd
from pathlib import Path

REQS_PER_SEC = 15


class GeneExpressionOnline:
    # stores basic class variables upon gene_expression_instance formation
    def __init__(self, gene_id):
        self.gene_id = gene_id
        self.req_count = 0
        self.last_req = 0

    # core function for getting expression data, calls other functions
    def get_expression(self):
        reference_gene_code_id, status_code = self.get_reference_gene_id()

        if status_code == 200:
            json_data, status_code = self.retrieve_expression_data(reference_gene_code_id)
            if status_code == 200:
                tissue_expression = {}
                for i in range(len(json_data["medianGeneExpression"])):
                    for k in json_data["medianGeneExpression"][i].keys():
                        tissue = json_data["medianGeneExpression"][i].get('tissueSiteDetailId')
                        expression_level = json_data["medianGeneExpression"][i].get("median")
                        # print(k+':', json_data["medianGeneExpression"][i][k])
                        tissue_expression[tissue] = expression_level
        else:
            tissue_expression = None
        return tissue_expression, status_code

    # as GTEx needs special (versioned?) gencode IDs this functions uses a special API endpoint to get the relevant gencodeID
    def get_reference_gene_id(self):

        # check if we need to rate limit ourselves
        if self.req_count >= REQS_PER_SEC:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        self.req_count += 1

        # try to pull the requested data and check if request was successful
        r = requests.get(
            "https://gtexportal.org/rest/v1/reference/gene?geneId=" + self.gene_id + "&gencodeVersion=v26&genomeBuild=GRCh38%2Fhg38&pageSize=250&format=json")

        if r.status_code == 200:
            json_data = r.json()
            response_list = json_data.get("gene")
            if response_list is not None:
                if len(response_list) > 1:
                    sys.stderr.write(
                        "There has been an ambiguity while looking for reference gencodeID!")
                if safe_get(response_list, 0):
                    return response_list[0].get("gencodeId"), 200
            else:
                return None, 494

        # if error "too many requests" retry after given time
        elif r.status_code == 429:
            if 'Retry-After' in r.headers:
                retry = r.headers['Retry-After']
                time.sleep(float(retry))
                self.retrieve_data()
        else:
            return None, r.status_code

    # function that takes in the reference gencodeID and retrieves expression data
    def retrieve_expression_data(self, reference_gene_code_id):

        # check if we need to rate limit ourselves
        if self.req_count >= REQS_PER_SEC:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        self.req_count += 1

        # try to pull the requested data and check if request was succesfull
        r = requests.get(
            "https://gtexportal.org/rest/v1/expression/medianGeneExpression?datasetId=gtex_v8&gencodeId=" + reference_gene_code_id + "&tissueSiteDetailId=Adipose_Subcutaneous%2CAdipose_Visceral_Omentum%2CAdrenal_Gland%2CArtery_Aorta%2CArtery_Coronary%2CArtery_Tibial%2CBladder%2CBrain_Amygdala%2CBrain_Anterior_cingulate_cortex_BA24%2CBrain_Caudate_basal_ganglia%2CBrain_Cerebellar_Hemisphere%2CBrain_Cerebellum%2CBrain_Cortex%2CBrain_Frontal_Cortex_BA9%2CBrain_Hippocampus%2CBrain_Hypothalamus%2CBrain_Nucleus_accumbens_basal_ganglia%2CBrain_Putamen_basal_ganglia%2CBrain_Spinal_cord_cervical_c-1%2CBrain_Substantia_nigra%2CBreast_Mammary_Tissue%2CCells_Cultured_fibroblasts%2CCells_EBV-transformed_lymphocytes%2CCells_Transformed_fibroblasts%2CCervix_Ectocervix%2CCervix_Endocervix%2CColon_Sigmoid%2CColon_Transverse%2CEsophagus_Gastroesophageal_Junction%2CEsophagus_Mucosa%2CEsophagus_Muscularis%2CFallopian_Tube%2CHeart_Atrial_Appendage%2CHeart_Left_Ventricle%2CKidney_Cortex%2CKidney_Medulla%2CLiver%2CLung%2CMinor_Salivary_Gland%2CMuscle_Skeletal%2CNerve_Tibial%2COvary%2CPancreas%2CPituitary%2CProstate%2CSkin_Not_Sun_Exposed_Suprapubic%2CSkin_Sun_Exposed_Lower_leg%2CSmall_Intestine_Terminal_Ileum%2CSpleen%2CStomach%2CTestis%2CThyroid%2CUterus%2CVagina%2CWhole_Blood&format=json")

        # if content is not empty return content
        if r.status_code == 200:
            return r.json(), 200
        elif r.status_code == 429:
            if 'Retry-After' in r.headers:
                retry = r.headers['Retry-After']
                time.sleep(float(retry))
                self.retrieve_expression_data(reference_gene_code_id)
        else:
            return None, r.status_code

class GeneExpression:
    def __init__(self, gene_id):
        self.expression_data = pd.read_csv(str(Path(__file__).parent) + "/assets/gtex_data.csv")
        self.ensemble_id = gene_id
        self.columns = ['Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)',
       'Adrenal Gland', 'Artery - Aorta', 'Artery - Coronary',
       'Artery - Tibial', 'Bladder', 'Brain - Amygdala',
       'Brain - Anterior cingulate cortex (BA24)',
       'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere',
       'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)',
       'Brain - Hippocampus', 'Brain - Hypothalamus',
       'Brain - Nucleus accumbens (basal ganglia)',
       'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)',
       'Brain - Substantia nigra', 'Breast - Mammary Tissue',
       'Cells - Cultured fibroblasts', 'Cells - EBV-transformed lymphocytes',
       'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid',
       'Colon - Transverse', 'Esophagus - Gastroesophageal Junction',
       'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube',
       'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex',
       'Kidney - Medulla', 'Liver', 'Lung', 'Minor Salivary Gland',
       'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary',
       'Prostate', 'Skin - Not Sun Exposed (Suprapubic)',
       'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum',
       'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina',
       'Whole Blood']

    def get_expression(self):
        gene_expression_data = self.expression_data.loc[self.expression_data.ensemble_id == self.ensemble_id]
        if len(gene_expression_data) == 0:
            return GeneExpressionOnline(self.ensemble_id).get_expression()
        else:
            return gene_expression_data[self.columns].to_dict(orient="records")[0], 200

# print(timeit.timeit(GeneExpression("ENSG00000198711").get_expression, number=50))
