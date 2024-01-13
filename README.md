
# PSPHunter: A Machine Learning Model to Predict Phase Separation Driving Residues

<div style="text-align: justify"> Dissecting the functions and the regulatory mechanisms of intracellular phase separation is fundamental to understanding transcriptional control, cell fate transition and disease development. However, the driving residues, which impact phase separation the most and therefore is the key for the functional study of protein phase separation, remain largely undisclosed. We developed PSPHunter, a machine learning method for predicting driving residues in phase-separating proteins. Validation through in vivo and in vitro methods, including FRAP and saturation measurements, confirms PSPHunter's accuracy. Applying PSPHunter, we demonstrate that truncating just 6 driving residues in SOX2 and GATA3 significantly disrupts their phase separation properties. Furthermore, PSPHunter identified nearly 80% of the phase-separating proteins associated with diseases. Remarkably, frequently mutated pathological residues (glycine and proline) tend to localize within driving residues, exerting a significant influence on phase separation. PSPHunter thus emerges as a crucial tool to uncover driving residues, facilitating insights into phase separation mechanisms governing transcriptional control, cell fate transitions, and disease development. </div>
--------------------------

![figs/overview.jpg](https://github.com/jsun9003/PSPHunter/blob/main/figs/overview.jpg)

## To generate the ensential features

#### Please have the following softwares installed first:
- psiblast for PSSM, https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- hhsuit for remote homology detection, https://github.com/soedinglab/hh-suite
- SPINE-D-local for Intrinsically disordered regions (IDRs) detection, http://zhouyq-lab.szbl.ac.cn/download/
- SNBRFinder for DNA and RNA binding prediction, http://ibi.hzau.edu.cn/SNBRFinder
- GPS5.0 for PTM prediction, https://gps.biocuckoo.cn/

## Additional Tutorial
- [hhsuit for HMM (Steinegger M. et al. BMC Bioinformatics, 2019)](https://github.com/soedinglab/hh-suite)
- [scikit-learn (Pedregosa, F. et al. Journal of Machine Learning Research, 2011)](https://scikit-learn.org/stable/getting_started.html)

## Datasets for Training
- the phase separation proteins used to contruct PSPHunter are in ./datasets
- 以训练好的模型存储在./train/里

## Generate the features
- 所有的feature产生的代码存放在scripts/featureExtraction, 包含序列特征和功能特征, 合并完成之后可用于模型的训练

## demo of predict phase separation probility
`cd Test`

`perl scripts/Standalone/predict_proteinProb.pl Test`


## Cite

Cite our paper by

```
@article {Sun2024multiscale,
	author = {Jun Sun, Jiale Qu and Cai Zhao},
	title = {Precise prediction of phase-separation key residues by machine learning reveals the dysregulated landscape impacted by pathogenic mutations},
	year={2024},
	publisher = {Nature Publishing Group},
	journal = {Under Submited}
}
```

![figs/graphAbstract.jpg](https://github.com/jsun9003/PSPHunter/blob/main/figs/graphAbstract.jpg)



# Availability

http://psphunter.stemcellding.org  

## Contact

Please contact o.sj@live.com or raise an issue in the github repo with any questions about installation or usage. 
