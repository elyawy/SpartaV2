import os
import re


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



def get_substitution_model(path):
    joined_path = os.path.join(path,"T1.raxml.log")

    res_dict = parse_raxmlNG_output(joined_path)

    x,y,z,m,n,k = res_dict["subAC"],res_dict["subAG"],res_dict["subAT"],res_dict["subCG"],res_dict["subCT"], res_dict["subGT"]
    x,y,z,m,n,k = float(x),float(y),float(z),float(m),float(n),float(k)
    pi_A,pi_C,pi_G,pi_T = res_dict["fA"],res_dict["fC"],res_dict["fG"],res_dict["fT"]
    pi_A,pi_C,pi_G,pi_T = float(pi_A),float(pi_C),float(pi_G),float(pi_T)

    a,b,c,d,e = (n/y)*(pi_G/pi_T), (z/y)*(pi_G/pi_T), (1.0/y)*(pi_G/pi_T), (x/y)*(pi_G/pi_C), (m/y)

    subtitution_model = {
        "freq": (pi_A, pi_C, pi_G, pi_T),
        "rates": (x,y,z,m,n,k),
        "inv_prop": res_dict["pInv"],
        "gamma_shape": res_dict["gamma"],
        "gamma_cats": res_dict["cats"],
        "submodel": "GTR"
    }
    return subtitution_model