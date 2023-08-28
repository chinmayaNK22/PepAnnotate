# PepAnnotate
Annotates multiple Peptide Spectrum Matches (PSMs) from proteomic database search result files

Basically, the `annotate_spectra` function under `pylab_aux` of pyteomics (version 4.5.6) package is used iteratively over PSMs.

### Inputs
Raw LC-MS/MS spectra in MGF/mzML format and database search result in pepxml format needs to be provided as inputs
```
# Usage

> python annotatespec.py [-h] [-r ...] -p [-p ...] -i [-i ...]

Plot Annotated Peptide Spectrum Matches (PSM)

positional arguments:
  -r          Path to folder where DDA derived raw files (MS/MS data)are stored in mzML or mgf format
  -p          PSMs from database search output in pep.xml format or Prosit predicted peptide spectral library in MSP format
  -i          Path to PSM file from Proteome Discoverer

options:
  -h, --help  show this help message and exit
```
