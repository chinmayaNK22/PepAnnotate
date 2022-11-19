from pyteomics import mgf, pepxml, mass, parser
import os
from urllib.request import urlopen, Request
import pylab
from pyteomics import pylab_aux as pa, usi
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from itertools import islice


inpath = "..\PRIDE\Mavium_Mavium_hominissuis_proteome_search_unassigned"

pepfile = "..\PRIDE\M_avium_MTBH37Rv_Ra_proteome_search_081821.pep.xml"

#infile = "M_avium_Mavium_hominissuis_proteome_search_081821_PSMs_filtered.txt"

def fetch_spectrum(infile):
    dicts = {}
    with mgf.read(infile) as spectra:
        for spectrum in spectra:
            spec_info = spectrum['params']['title']
            try:
                file_path = spec_info.split(',')[0].split('File:')[1].strip('"')
                raw_file = os.path.split(file_path)[-1]
                scan = spec_info.split(',')[-1].split('scan=')[1].strip('"')

            except:
                file_path = spec_info.split(';')[0].lstrip('File: ').strip('"')
                raw_file = os.path.split(file_path)[-1]
                scan = spec_info.split(';')[-1].lstrip('scans: ').strip('"')

            if raw_file + '@' + scan not in dicts:
                dicts[raw_file + '@' + scan] = [spectrum]
            else:
                dicts[raw_file + '@' + scan].append(spectrum)
                
    rawfile_format = [k.split('@')[0].split('.')[-1] for k, v in dicts.items()]
    return dicts, rawfile_format[0]

def read_mgf_files(inpath):
    infiles = []
    for files in os.listdir(os.path.join(inpath)):
        if os.path.isfile(os.path.join(inpath, files)):
            if files.split('.')[-1] == 'mgf':
                infiles.append(os.path.join(inpath, files))
    return infiles

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

def get_modX_pep(mod_pep, mods):
    if len(mods) == 1:
        stripped_seq = mod_pep.split('.')[1]
        mod_pos = int(mods[0]['position']) - 1
        mod_aa = str(round(mods[0]['mass']))+stripped_seq[mod_pos].upper()
        if parser.is_modX(mod_aa) == True:
            modx_seq = stripped_seq[0:mod_pos] + mod_aa + stripped_seq[mod_pos+1::]
            return modx_seq
    else:
        print (mod_pep, mods)


def read_psm_txt_file(infile):
    gssp = {}
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            pep = split_i[4].split('.')[1].upper()
            pro = split_i[8].strip('"')
            category = ""
            if len(category) != 0:
                if pep not in gssp:
                    gssp[pep] = [category + '@' + pro]
                else:
                    gssp[pep].append(category + '@' + pro)
            else:
                if pep not in gssp:
                    gssp[pep] = [pro]
                else:
                    gssp[pep].append(pro)

    return gssp


def plot_spectra(spec_file, scan, all_spectra, psm_type, peps, mod_pep, rt, score_type, score, file_format):

    rawfile = ""
    if spec_file.split('.')[-1] != file_format:
        rawfile = spec_file.rstrip(spec_file.split('.')[-1]) + file_format
    else:
        rawfile = spec_file

    try:
        if rawfile + '@' + scan in all_spectra:
            spectrum = all_spectra[rawfile + '@' + scan][0]
            #pylab.figure()
            #pylab_aux.plot_spectrum(spectrum, title="Experimental spectrum " + spectrum['params']['title'])

            fig, ax = plt.subplots(figsize=(12, 6))
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            if len(psm_type) != 0:
                fig_name = psm_type + '_' + peps['search_hit'][0]['peptide'] + '_' + peps['spectrum']
                header = psm_type + ', ' + 'Sequence: ' + mod_pep + ', RT (min): ' + rt + ', ' + score_type + ': ' + score
                pa.annotate_spectrum(spectrum, peps['search_hit'][0]['peptide'], title= header, maxcharge=peps['assumed_charge'], ion_types='aby')
                #plt.show()     
                plt.savefig(fig_name + '.pdf', bbox_inches="tight", dpi=300)
                plt.close()
                #plt.savefig(fig_name + '.png', bbox_inches="tight", dpi=300)
            else:
                fig_name = peps['search_hit'][0]['proteins'][0]['protein'] + '_' + peps['search_hit'][0]['peptide'] + '_' + peps['spectrum']
                header = 'Sequence: ' + mod_pep + ', RT (min): ' + rt + ', ' + score_type + ': ' + score
                pa.annotate_spectrum(spectrum, peps['search_hit'][0]['peptide'], title= header, maxcharge=peps['assumed_charge'], ion_types='aby')
                #plt.show()     
                plt.savefig(fig_name + '.pdf', bbox_inches="tight", dpi=300)
                plt.close()
                #plt.savefig(fig_name + '.png', bbox_inches="tight", dpi=300)

    except:
        print ('ERROR: Could not annotate peptide ', mod_pep)
   
def annotate_peptide_spectra(mgf_path, pepfile, psm_file):
    
    psms = fetch_pepxml_data(pepfile)

    sorted_psms = {}
    if len(psm_file) != 0:
        sorted_psms = read_psm_txt_file(psm_file)
    
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

    for infile in read_mgf_files(inpath):
        if os.path.split(infile)[-1] in psms:
            try:
                all_spectra, file_format = fetch_spectrum(infile)
                
                print ('Stored ', len(all_spectra), ' MS/MS spectra from ', os.path.split(infile)[-1])
                
                for peps in high_score_psms:
                    rt = str(round(float(peps['retention_time_sec']/60),2))
                    z = peps['assumed_charge']
                    mod_pep = peps['search_hit'][0]['modified_peptide']
                    score_type = next(iter(peps['search_hit'][0]['search_score']))
                    score = str(round(peps['search_hit'][0]['search_score'][score_type],2))
                    mods = peps['search_hit'][0]['modifications']
                    peptide = ""
                    if len(mods) != 0:
                        peptide = get_modX_pep(mod_pep, mods)
                    else:
                        peptide = peps['search_hit'][0]['peptide']
                        
                    spec_file = peps['spectrum'].split('.')[0] + '.mgf'
                    scan = str(peps['start_scan'])
        
                    if len(psm_file) != 0:
                        try:
                            for sorted_peps in sorted_psms[peps['search_hit'][0]['peptide']]:
                                psm_type = ""
                                prot = ""
                                if '@' in sorted_peps:
                                    psm_type = sorted_peps.split('@')[0]
                                    prot = sorted_peps.split('@')[1]
                                else:
                                    psm_type = sorted_peps

                                plot_spectra(spec_file, scan, all_spectra, psm_type, peps, mod_pep, rt, score_type, score, file_format)
          
                        except:
                            pass
                            #print ('Sequence ', peps['search_hit'][0]['peptide'], ' is not there in ', os.path.split(psm_file)[-1])

                    else:
                        plot_spectra(spec_file, scan, all_spectra, "", peps, mod_pep, rt, score_type, score, file_format)

                        
            except:
                print ('Could not process annotate spectra from the file ', os.path.split(infile)[-1])

                
annotate_peptide_spectra(inpath, pepfile, "")
        
        
