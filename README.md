
# PSPHunter: A Machine Learning Model to Predict Phase Separation Driving Residues

Dissecting the functions and the mechanisms of intracellular phase separation regulation is fundamental to understanding transcriptional control, cell fate transition and disease development. However, the limited knowledge of driving residues, which most impact phase separation, limits the functional study of protein phase separation, particularly in terms of how abnormal phase separation leads to disease. Here, we developed PSPHunter, a machine learning method that can precisely predict driving residues of phase-separating proteins. Through systematically analyzing and experimentally verifying, we confirmed that PSPHunter can predict known driving residues of typical phase-separating proteins and new driving residues of SOX2 and GATA3. Notably, truncation of only 6 driving residues could significantly disrupt their phase separation abilities. Furthermore, we find that nearly 80% of phase-separating proteins were highly disease-related. Most frequent pathogenic mutations of phase-separating proteins, glycine and proline, are preferentially located at the driving residues and significantly impact phase separation than other mutations. Altogether, we illustrated that PSPHunter could provide an analytic framework for dissecting the relationship among amino acid residues, phase separation and diseases.
--------------------------

![figs/overview.jpg](https://github.com/jsun9003/scCARE-seq/blob/main/figs/overview.jpg)

# Preprocessing of scCARE-seq datasets

#### Please have the following softwares installed first:
- bowtie2, http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
- samtools, http://www.htslib.org/
   samtools version >= 1.3.1 is required.
- Trim_galore, https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
- STAR, https://github.com/alexdobin/STAR
- Higashi, https://github.com/ma-compbio/Higashi
- sci-CAR, https://github.com/JunyueC/sci-CAR_analysis

- Optional: FastQC, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Additional Tutorial
- [Higashi-analysis for HiC (Zhang et al.Nature biotechnology, 2022)](https://github.com/ma-compbio/Higashi)
- [sci-CAR_analysis for RNA (Cao et al.Science, 2018)](https://github.com/JunyueC/sci-CAR_analysis)

# Analysis of scCARE-seq datasets include the following steps:
## 1. Split scCARE-seq data into HiC partion and RNA partion

## 2. Single cell HiC analysis for the HiC partion

## 3. Single cell RNA-seq analysis for the RNA partion

# Cite

Cite our paper by

```
@article {Sun2022multiscale,
	author = {Jun Sun, Jiale Qu and Cai Zhao},
	title = {PSPHunter: A Machine Learning Model to Predict Phase Separation Driving Residues Enriching Pathogenic Mutations of Glycine and Proline},
	year={2022},
	publisher = {Nature Publishing Group},
	journal = {Under Submited}
}
```

![figs/graphAbstract.jpg](https://github.com/jsun9003/scCARE-seq/blob/main/figs/graphAbstract.jpg)



# Contact

Please contact o.sj@live.com or raise an issue in the github repo with any questions about installation or usage. 
