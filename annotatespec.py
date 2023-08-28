from pyteomics import mgf, pepxml, mass, parser, mzml
import os
from pyteomics import pylab_aux as pa, usi
from itertools import islice
from src import parse_msp
from src import parse_psm_scan
from src import plot
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='''Plot Annotated Peptide Spectrum Matches (PSM)''')

parser.add_argument('inpath', metavar='-r', type=str, nargs='*', help='Path to folder where DDA derived raw files (MS/MS data)are stored in mzML or mgf format')

parser.add_argument('pepfile', metavar='-p', type=str, nargs='+', help='PSMs from database search output in pep.xml format or Prosit predicted peptide spectral library in MSP format')

parser.add_argument('infile', metavar='-i', type=str, nargs='+', help='Path to PSM file from Proteome Discoverer')

#parser.add_argument('tol', metavar='-t', type=str, nargs='+', help="Set the fragment tolerance in 'ppm' or 'Da' for matching experimental and theoretical MS/MS peaks")

args = parser.parse_args()

#inpath = "E:\\NTMs_Raw_Data\M_abscessus\In-house\Raw_Data"

#pepfile = "E:\\NTMs_Raw_Data\MSMS_Cluster\SSPP_Prosit_MSPs\M_abscessus_SSPPs_SpecLib.msp"

#infile = ["E:\\NTMs_Raw_Data\MSMS_Cluster\PSM_RawScan_Mapped\M_abscessus_In_house_specific_rawfile_cluster_psm_rawinfo_mapped.txt"]

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

def store_rawfiles(inpath):
    fformat = {'mgf':True, 'mzML':True, 'RAW':False, 'wiff':False}
    infiles = {}
    for files in os.listdir(os.path.join(inpath)):
        if os.path.isfile(os.path.join(inpath, files)):
            if files.split('.')[-1] in fformat:
                if fformat[files.split('.')[-1]] == True:
                    if files not in infiles:
                        infiles[files] = [os.path.join(inpath, files)]
                        
                else:
                    raise ValueError(f'ERROR: The input file in not MGF or mzML format. Convert the raw data into MGF or mzML format using ProteoWizard.')
                
    return infiles


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

def fetch_precursor_info(scan_num, scan):
    if 'precursorList' in scan:
        precursor_info = scan['precursorList']['precursor'][0]['selectedIonList']
        mz = precursor_info['selectedIon'][0]['selected ion m/z']
        try:
            z = precursor_info['selectedIon'][0]['charge state']
        except:
            #print (precursor_info['selectedIon'][0])
            z = "-"
        try:
            intensity = precursor_info['selectedIon'][0]['peak intensity']
        except:
            intensity = '0'
            print (f'WARNING: Count not find precursor intensity for {mz} of {scan_num}. Imputing it with value 0')
            #scan_n = precursor_info['spectrumRef']

        rt = ''
        base_mz = ''
        scan_info = scan['scanList']['scan']
        for scan_details in scan_info:
            rt = scan_details['scan start time']
            base_mz = scan['base peak m/z']
            
        return mz, rt, z, intensity
    
    else:
        #print (f'{scan_num} is a precursor scan')
        rt =''
        base_mz = ''
        scan_info = scan['scanList']['scan']
        for scan_details in scan_info:
            rt = scan_details['scan start time']
            base_mz = scan['base peak m/z']

        return 'Precursor_scan', rt, '-', '0'

def read_psm_txt_file(infile):
    peptides = {}
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            pep = split_i[4].split('.')[1].upper()
            pro = split_i[8].strip('"')
            category = ""
            if len(category) != 0:
                if pep not in peptides:
                    peptides[pep] = [category + '@' + pro]
                else:
                    peptides[pep].append(category + '@' + pro)
            else:
                if pep not in peptides:
                    peptides[pep] = [pro]
                else:
                    peptides[pep].append(pro)

    return peptides

def L2_norm(intensity_array):
    '''Performs L2 normalization of an intensity array.
    Refered the code from https://github.com/wilhelm-lab/spectrum_fundamentals/blob/development/spectrum_fundamentals/metrics/similarity.py'''

    #L2 normalization
    norm_intensity = np.linalg.norm(intensity_array, 2)

    # Normalize the vector to unit length
    normalized_vector = intensity_array / norm_intensity

    return normalized_vector

def map_peptide_to_spectra(msp_file, infiles):

    psms = parse_psm_scan.psm_rawinfo(infiles)

    predicted_pep_spec = parse_msp.extract_msp(msp_file)

    mapped_psms = 0
    pep_raw_scan_mapped = {}
    for predicted_pep in predicted_pep_spec:
        peptide = predicted_pep[0]
        z = predicted_pep[1]
        attrbs = predicted_pep[2]
        pred_mzs = predicted_pep[3]
        pred_intensities = predicted_pep[5]
        pred_ions = predicted_pep[4]

        nr_raw_scans = {}
        if peptide in psms:
            mapped_psms += 1
            for rawfile_scan in psms[peptide]:
                if z == rawfile_scan.split('@')[2]:
                    exp_rawfile = rawfile_scan.split('@')[0] + '.mzML'
                    exp_scan = rawfile_scan.split('@')[1]
                    nr_raw_scans[exp_rawfile + '@' + exp_scan + '@' + peptide + '@' + z] = 0

        for f, s in nr_raw_scans.items():
            rawf = f.split('@')[0]
            raw_scan = f.split('@')[1]
            pep = f.split('@')[2]
            charge = f.split('@')[3]
            
            if rawf not in pep_raw_scan_mapped:
                pep_raw_scan_mapped[rawf] = [[raw_scan + '@' + pep + '@' + charge] + [pred_mzs] + [pred_intensities] + [pred_ions]]
            else:
                pep_raw_scan_mapped[rawf].append([raw_scan + '@' + peptide + '@' + z] + [pred_mzs] + [pred_intensities] + [pred_ions])

    return pep_raw_scan_mapped
  
def annotate_peptide_spectra(rawfile_path, pepfile, psm_file):

    rawfiles = store_rawfiles(rawfile_path)

    #def psm_msp_map():
    if pepfile.split('.')[-1] == 'zip' or pepfile.split('.')[-1] == 'msp':
        mapped_peptides = map_peptide_to_spectra(pepfile, psm_file)

        for rawfile, info in mapped_peptides.items():
            try:
                r_file = rawfiles[rawfile][0]
            except:
                raise ValueError(f'ERROR: Raw file {rawfile} was not found in the input folder {rawfile_path}.')
            
            rawinfo = mzml.read(r_file, read_schema=True)
            scans = {rinfo['id'].split(' ')[-1].strip('scan='):rinfo for rinfo in rawinfo}

            print (f'INFO: Stored spectra of {len(info)} peptide spectrum matches.')
                
            scan_count = 0
            for scan_pep in info:
                
                exp_scan_num = scan_pep[0].split('@')[0]
                sequence = scan_pep[0].split('@')[1]
                charge = scan_pep[0].split('@')[2]
                
                if exp_scan_num in scans:
                    precursor, rt, z, mz_intensity = fetch_precursor_info(exp_scan_num, scans[exp_scan_num])
                    print (precursor, rt, z, mz_intensity)
                    exp_spectra = scans[exp_scan_num]
                    spec_mz = exp_spectra['m/z array']
                    spec_intensity = exp_spectra['intensity array']
                    pred_mz = scan_pep[1]
                    pred_intensity = scan_pep[2]
                    pred_ions = scan_pep[3]

                    #print (r_file, exp_scan_num, sequence, charge)

                    # Normalize the intensity array
                    norm_spec_intensity = L2_norm(spec_intensity)
                    norm_pred_intensity = L2_norm(pred_intensity)
                    
                    plot.plot_mirror_plot(r_file, exp_scan_num, sequence, precursor, rt, mz_intensity, charge, spec_mz, norm_spec_intensity, pred_mz, norm_pred_intensity, pred_ions)

    elif 'pep.xml' in pepfile:
        high_score_psms, psms = manipulate_pepxml(pepfile)

        sorted_psms = {}
        if len(psm_file) != 0:
            sorted_psms = read_psm_txt_file(psm_file)
        
        for infile in rawfiles:
            if os.path.split(infile)[-1] in psms:
                try:
                    all_spectra, file_format = fetch_spectrum(infile)
                    
                    print (f'INFO: Stored {len(all_spectra)} MS/MS spectra from {os.path.split(infile)[-1]}')
                    
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

                                    plot.plot_spectra(spec_file, scan, all_spectra, psm_type, peps, mod_pep, rt, score_type, score, file_format)
              
                            except:
                                pass
                                #print ('Sequence ', peps['search_hit'][0]['peptide'], ' is not there in ', os.path.split(psm_file)[-1])

                        else:
                            plot.plot_spectra(spec_file, scan, all_spectra, "", peps, mod_pep, rt, score_type, score, file_format)

                            
                except:
                    print (f'INFO: Could not process annotate spectra from the file {os.path.split(infile)[-1]}')
                    

if __name__== "__main__":

    if args.pepfile[0].split('.')[-1] == 'pep.xml':
        if len(args.infile) == 0:
            print (f'WARNING: No input PSM file was found.')
            annotate_peptide_spectra(args.inpath[0], args.pepfile[0], "")
        else:
            annotate_peptide_spectra(args.inpath[0], args.pepfile[0], args.infile)
        
    elif args.pepfile[0].split('.')[-1] == 'zip':
        annotate_peptide_spectra(args.inpath[0], args.pepfile[0], args.infile)

    elif args.pepfile[0].split('.')[-1] == 'msp':
        annotate_peptide_spectra(args.inpath[0], args.pepfile[0], args.infile)
    
        
        
