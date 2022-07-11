import requests
from tenacity import retry, stop_after_attempt, wait_exponential, wait_random

@retry(reraise=True,
       stop=stop_after_attempt(5),
       wait=wait_exponential(multiplier=1.3, min=0.1, max=3))
def vep_variant_recoder_request(variant):
    """General function for handling API communication. If there is some error with the returned data, the request
    will be retried a couple of times. Return obtained data as a dict.
    """
    url = f"http://grch37.rest.ensembl.org/variant_recoder/human/{variant}?content-type=application/json"

    # try to pull the requested data and check if request was successful
    r = requests.get(url, headers={"content-type": "application/json"})

    if r.ok:
        try:
            return r.json()[0], 200
        except IndexError:
            print("Error decoding JSON received from variant_recoder API!")
            return None, 200
    else:
        # Return None if some sort of different error occurs.
        print(f"VEP VARIANT RECODER ERROR '{r.status_code}: {r.reason}' occured for {variant}. Retrying...")
        raise IOError("There has been an issue with a variant.")

def convert_variant(variant):
    try:
        response_decoded, status_code = vep_variant_recoder_request(variant)
    except IOError:
        response_decoded = None
        status_code = 400

    if status_code == 200:
        try:
            hgvsc = response_decoded.get(list(response_decoded)[0]).get("hgvsc")
            for i, _item in enumerate(hgvsc):
                if ("ENST" not in _item) and ("NM" not in _item) or (":c." not in _item):
                    hgvsc.pop(i)
            return hgvsc
        except (IndexError, AttributeError):
            return []
