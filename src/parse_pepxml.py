from pyteomics import pepxml

def fetch_pepxml_data(pepfile):
    psm_dicts = {}
    with pepxml.read(pepfile) as psms:
        for psm in psms:
            spec_file = psm['spectrum'].split('.')[0] + '.mgf'
            if spec_file not in psm_dicts:
                psm_dicts[spec_file] = [psm]
            else:
                psm_dicts[spec_file].append(psm)
                
    return psm_dicts

def manipulate_pepxml(pepfile):
    psms = fetch_pepxml_data(pepfile)
    
    peptide_scores = {}
    psm_info = {}
    for file, psm in psms.items():
        for peps in psm:
            rt = str(float(peps['retention_time_sec']/60))
            z = peps['assumed_charge']
            mod_pep = peps['search_hit'][0]['modified_peptide']
            score_type = next(iter(peps['search_hit'][0]['search_score']))
            score = peps['search_hit'][0]['search_score'][score_type]
            if mod_pep not in peptide_scores:
                peptide_scores[mod_pep] = [score]
            else:
                peptide_scores[mod_pep].append(score)

            if mod_pep not in psm_info:
                psm_info[mod_pep] = [peps]
            else:
                psm_info[mod_pep].append(peps)

    high_score_psms = []
    for k, values in peptide_scores.items():
        high_score = sorted(values)[-1]
        idx_n = ''.join([str(idx) for idx, v in enumerate(values) if high_score == v])
        high_score_psms.append(psm_info[k][int(idx_n)])

    print ('There are ', len(high_score_psms), ' PSMs stored for spectrum annotation.')

    return high_score_psms, psms
