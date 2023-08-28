import matplotlib.pyplot as plt
from pyteomics import mass
from pyteomics import pylab_aux as pa, usi
import numpy as np
import os
import spectrum_utils.spectrum as sus
import spectrum_utils.iplot as sup

#Examples
peptide = 'DLTDYLoxMK'  # oxidized methionine
aa_mass = mass.std_aa_mass.copy()
aa_mass['ox'] = 15.9949  # define the mass of the label

##usi_top = "mzspec:PXD004732:01650b_BC2-TUM_first_pool_53_01_01-3xHCD-1h-R2:scan:41840"
####usi_bottom = 'mzspec:MSV000080679:j11962_C1orf144:scan:10671'
##
##print (usi_top)
##spectrum_top = usi.proxi(usi_top, 'massive')
##spectrum_bottom = usi.proxi(usi_bottom, 'massive')
##
#print (spectrum_top)

#spectrum_top1 = sus.MsmsSpectrum.from_usi(usi_top)
##spectrum_top1.annotate_proforma(peptide, 0.5, "Da", ion_types="aby")
##

def annotate_ion_type(annotation, ion_types="aby"):
    if annotation.ion_type[0] in ion_types:
        if abs(annotation.isotope) == 1:
            iso = "+i" if annotation.isotope > 0 else "-i"
        elif annotation.isotope != 0:
            iso = f"{annotation.isotope:+}i"
        else:
            iso = ""
        nl = {"-NH3": "*", "-H2O": "o"}.get(annotation.neutral_loss, "")
        return f"{annotation.ion_type}{iso}{'+' * annotation.charge}{nl}"
    else:
        return ""   

def spec_utils_compatible(rawfile:str, scan_num:str, precursor:str, rt:str, mz_intensity:str, charge:str, mz: np.ndarray, intensity: np.ndarray):
    '''Assemble a dictionary which is compatible as an input to the pyteomics mirror plot with MS/MS spectra, its corresponding intensity array with other details'''

    spec_utils_input = {'attributes':{"precursor_mz":0,"precursor_charge":0}, 'm/z array':0,'intensity array':0,'status':str,'usi':str}

    spec_utils_input['attributes']["precursor_mz"] = precursor
    spec_utils_input['attributes']["precursor_charge"] = charge
    
    spec_utils_input['m/z array'] = mz
    
    spec_utils_input['intensity array'] = intensity
    
    spec_utils_input['status'] = 'READABLE'
    
    rawinfo = 'mzspec:' + 'In-house' + ':ccms_peak' + rawfile + ':scan:' + scan_num
    
    spec_utils_input['usi'] = rawinfo

    
    #spectrum_top1 = sus.MsmsSpectrum.from_usi(spec_utils_input, precursor_mz = precursor, precursor_charge = charge)
    #print (spectrum_top1)
    
    return spec_utils_input

def plot_mirror(rawfile, exp_scan_num, spectrum_top: dict, spectrum_bottom: dict, peptide: str, charge: int):

    '''Plot two mass spectra of a same peptide as peptide fragments annotated mirror plot''' 
    
    figure_name = ':'.join([peptide, str(charge), os.path.split(rawfile)[-1], exp_scan_num])
                           
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.spines["right"].set_visible(True)
    ax.spines["top"].set_visible(True)
    pa.mirror(spectrum_top, spectrum_bottom, peptide=peptide, precursor_charge=charge, aa_mass=aa_mass, annot_fmt=annotate_ion_type, ion_types='aby', ftol=0.05, grid = False, xlabel = 'm/z', ylabel = '% Rel. intensity', title = figure_name,
             scaling=None, remove_precursor_peak=True, backend='spectrum_utils', ax = ax)

    #chart = sup.mirror(spectrum_top, spectrum_bottom)
    
    figure_name = '-'.join([peptide, str(charge), os.path.split(rawfile)[-1], exp_scan_num]) + '.pdf'
    #plt.savefig(figure_name, dpi=300, bbox_inches="tight", transparent=True)
    print (f'INFO: Plotted mirror plot for the peptide {peptide} with {charge} into {figure_name}.')
    plt.show()
    #plt.close()
    

def plot_mirror_plot(rawfile, exp_scan_num, sequence, precursor, rt, mz_intensity, charge, spec_mz, spec_intensity, pred_mz, pred_intensity, pred_ions):

    '''Split the peptide annotated spectra from two sources as top and bottom spectra for the mirror plot generation'''
    experiment_spectra = spec_utils_compatible(rawfile, exp_scan_num, precursor, rt, mz_intensity, charge, spec_mz, spec_intensity)
    predicted_spectra = spec_utils_compatible('Prosit', 'Predicted', precursor, rt, mz_intensity, charge, pred_mz, pred_intensity)
    
    plot_mirror(rawfile, exp_scan_num, experiment_spectra, predicted_spectra, sequence, int(charge))


def plot_spectra(spec_file, scan, all_spectra, psm_type, peps, mod_pep, rt, score_type, score, file_format):

    '''Plot a peptide spectrum match with its fragment annotation using 'annotate_spectrum' function of pyteomics package'''
    
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
