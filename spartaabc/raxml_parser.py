import os
import re
from pathlib import Path

from msasim import sailfish as sf

MODEL_MAPPER = {
    # Protein models
    "WAG": sf.MODEL_CODES.WAG,
    "LG": sf.MODEL_CODES.LG,
    "cpREV": sf.MODEL_CODES.CPREV45,
    "Dayhoff": sf.MODEL_CODES.DAYHOFF,
    "HIVb": sf.MODEL_CODES.HIVB,
    "HIVw": sf.MODEL_CODES.HIVW,
    "JTT": sf.MODEL_CODES.JONES,
    "mtREV": sf.MODEL_CODES.MTREV24,
    # Nucleotide models
    "JC": sf.MODEL_CODES.NUCJC,
    "HKY": sf.MODEL_CODES.HKY,
    "GTR": sf.MODEL_CODES.GTR,


}

def parse_raxmlNG_output(res_filepath):

    try:
        with open(res_filepath) as fpr:
            content = fpr.read()
        res_dict = parse_raxmlNG_content(content)
    except:
        print("Error with:", res_filepath)
        return

    return res_dict


def parse_raxmlNG_content(content):
    """
    :return: dictionary with the attributes - string typed. if parameter was not estimated, empty string
    """
    res_dict = dict.fromkeys(["ll", "pInv", "gamma", "cats",
                              "fA", "fC", "fG", "fT",
                              "subAC", "subAG", "subAT", "subCG", "subCT", "subGT",
                              "time"], "")

    # likelihood
    ll_re = re.search("Final LogLikelihood:\s+(.*)", content)
    if ll_re:
        res_dict["ll"] = ll_re.group(1).strip()
    elif re.search("BL opt converged to a worse likelihood score by", content) or re.search("failed", content):
        ll_ini = re.search("initial LogLikelihood:\s+(.*)", content)
        if ll_ini:
            res_dict["ll"] = ll_ini.group(1).strip()
    else:
        res_dict["ll"] = 'unknown raxml-ng error, check "parse_raxmlNG_content" function'


    # gamma (alpha parameter) and proportion of invariant sites
    gamma_regex = re.search("alpha:\s+(\d+\.?\d*)\s+", content)
    pinv_regex = re.search("P-inv.*:\s+(\d+\.?\d*)", content)
    cats_regex = re.search("\(\d+ cats,", content)
    if gamma_regex:
        res_dict['gamma'] = gamma_regex.group(1).strip()
    if pinv_regex:
        res_dict['pInv'] = pinv_regex.group(1).strip()
    if cats_regex:
        res_dict['cats'] = cats_regex.group(0)[1]

    # Nucleotides frequencies
    nucs_freq = re.search("Base frequencies.*?:\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)", content)
    if nucs_freq:
        for i,nuc in enumerate("ACGT"):
            res_dict["f" + nuc] = nucs_freq.group(i+1).strip()

    # substitution frequencies
    subs_freq = re.search("Substitution rates.*:\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)", content)
    if subs_freq:
        for i,nuc_pair in enumerate(["AC", "AG", "AT", "CG", "CT", "GT"]):  # todo: make sure order
            res_dict["sub" + nuc_pair] = subs_freq.group(i+1).strip()

    # Elapsed time of raxml-ng optimization
    rtime = re.search("Elapsed time:\s+(\d+\.?\d*)\s+seconds", content)
    if rtime:
        res_dict["time"] = rtime.group(1).strip()
    else:
        res_dict["time"] = 'no ll opt_no time'

    return res_dict




def parse_raxml_bestModel(model_path: Path):
    try:
        model_file = next(model_path.glob("*.bestModel"))
    except StopIteration as e:
        raise Exception(f"No substitution model provided in {model_path}, please specify a model using the RaxML-ng bestModel format.")

        # Initialize results dictionary
    results = {
        'submodel': None,
        'empirical_frequencies': False,
        'gamma_cats': 1,
        'gamma_shape': 1.0,
        'partition': None
    }
    
    # Split into model and partition if partition exists
    parts = model_file.read_text().split(',')
    model_part = parts[0].strip()
    
    # Parse partition information if present
    if len(parts) > 1:
        partition_part = parts[1].strip()
        results['partition'] = partition_part
    
    # Split model into components
    model_components = model_part.split('+')
    
    # First component is always the substitution model
    left_bracket_index = model_components[0].find("{")
    if left_bracket_index == -1:
        left_bracket_index = None

    results['submodel'] = MODEL_MAPPER.get(model_components[0][:left_bracket_index], -1)
    if (results['submodel']) == -1:
        raise RuntimeError("The requested model has not been implemented :(")
    
    if results['submodel'] == sf.MODEL_CODES.GTR:
        left_bracket_index = model_components[0].find("{")
        right_bracket_index = model_components[0].find("}")
        rates = model_components[0][left_bracket_index+1:right_bracket_index]
        rates = rates.split("/")
        rates = [float(x) for x in rates]
        results["params"] = rates


    # Parse remaining components
    for component in model_components[1:]:
        # Check for empirical frequencies
        left_bracket_index = component.find("{")
        right_bracket_index = component.find("}")

        if component.startswith('F'):
            rates = []
            if left_bracket_index != -1:
                frequencies = (component[left_bracket_index+1:right_bracket_index])  # Remove 'G' and convert to int
                frequencies = frequencies.split("/")
                frequencies = [float(x) for x in frequencies]

                results['params'] = frequencies + results['params']
        
        # Check for gamma categories and alpha
        elif component.startswith('G'):
            # Extract number of categories
            # gamma_info = component.split('m')
            if "m" in component:
                shift = -1

            results['gamma_cats'] = int(component[1:left_bracket_index+shift])  # Remove 'G' and convert to int
            
            # Extract alpha if present
            if left_bracket_index != -1:
                alpha_str = component[left_bracket_index:].strip('{}')
                results['gamma_shape'] = float(alpha_str)

        elif component.startswith('I'):
            left_bracket_index = component.find("{")
            right_bracket_index = component.find("}")
            if left_bracket_index != -1:
                alpha_str = component[left_bracket_index:].strip('{}')
                results['invariant_sites'] = float(alpha_str)


    
    return results
    



def get_substitution_model(path):
    return parse_raxml_bestModel(path)
    # joined_path = os.path.join(path,"T1.raxml.log")

    # res_dict = parse_raxmlNG_output(joined_path)

    # x,y,z,m,n,k = res_dict["subAC"],res_dict["subAG"],res_dict["subAT"],res_dict["subCG"],res_dict["subCT"], res_dict["subGT"]
    # x,y,z,m,n,k = float(x),float(y),float(z),float(m),float(n),float(k)
    # pi_A,pi_C,pi_G,pi_T = res_dict["fA"],res_dict["fC"],res_dict["fG"],res_dict["fT"]
    # pi_A,pi_C,pi_G,pi_T = float(pi_A),float(pi_C),float(pi_G),float(pi_T)

    # a,b,c,d,e = (n/y)*(pi_G/pi_T), (z/y)*(pi_G/pi_T), (1.0/y)*(pi_G/pi_T), (x/y)*(pi_G/pi_C), (m/y)

    # subtitution_model = {
    #     "freq": (pi_A, pi_C, pi_G, pi_T),
    #     "rates": (x,y,z,m,n,k),
    #     "inv_prop": res_dict["pInv"],
    #     "gamma_shape": res_dict["gamma"],
    #     "gamma_cats": res_dict["cats"],
    #     "submodel": "GTR"
    # }
    # return subtitution_model