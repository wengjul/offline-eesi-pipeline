# Offline Extractive Electrospray Mass Spectrometry data analysis pipeline

Zenodo DOI: https://doi.org/10.5281/zenodo.20542289


Code repository accompanying:


Julian Weng et al. (2026)  
*A Decadal-Scale Perspective on PM10 Composition and its Variability Drivers at the Alpine High-Altitude Research Station Jungfraujoch*  
Atmospheric Chemistry and Physics


## Authors

Julian Weng, Yufang Hao, Patrik Winiger, Imad El Haddad, and Kaspar R. Daellenbach

PSI Center for Energy and Environmental Sciences  
Paul Scherrer Institute (PSI)  
5232 Villigen PSI, Switzerland




## Overview

This repository contains Python workflows for processing and analyzing
offline Extractive Electrospray Mass Spectrometry (EESI-MS) measurements, as used in the manuscript.
The following analysis part is covered: from fitted HR EESI data (Tofware output) until aggregated scaled relative ion intensity time series, including QC steps. 
For details, refer to the SI of the manuscript.



The workflow includes:

1. Raw data import and preprocessing (including filtering steps)
2. Assignment of sample and waterblank periods
3. Normalization
4. Averaging and water blank subtraction
5. 'Semi-Quantification'
6. field blank subtraction
7. Additional steps, including QC and ion filtering (based on different criteria)




---



### Directories

| Directory | Description |
|------------|------------|
| `data example/` | Example input data |
| `output/` | Example processed data products |
| `workflow_overview_slides/` | Overview slides for the pipeline (as orientation, not reviewed) |
| `auxiliary_codes/` | Auxiliary codes, which are not used in the pipeline but essential for data quality checks and further processing (ion filtering) |

---

## Requirements

The code was developed using:

- Python 3.10.9
- NumPy
- Pandas
- Plotly



Create the environment:

```bash
conda env create -f environment-PM10-composition-JFJ.yaml
conda activate environment-PM10-composition-JFJ
```

---

## Data

### Input data

The repository uses:

- Offline Extractive Electrospray Mass Spectrometry (EESI-MS) measurements, after export from Tofware (after HR peak fitting), including files for MS data, ion formula and mass list, measurement time series and run numbers
- Ancillary offline autosampler data (timings)
- Additional sample information data (including e.g. exact extraction volumns)


Example data are inlcuded into the repository.
Full data are available upon justified request.

---


## Citation

If you use this code, please cite:

Julian Weng et al. (2026)

and the repository DOI, linked to Zenodo:

DOI: https://doi.org/10.5281/zenodo.20542289

---

## License

MIT License

See `LICENSE` for details.

---

## Contact

Julian Weng
Paul Scherrer Institute
julian.weng@psi.ch









