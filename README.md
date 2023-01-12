# A Network medicine approach to study comorbidities in heart failure with preserved ejection fraction

## Abstract

**Background**:  
Comorbidities are expected to impact pathophysiology of heart failure (HF) with preserved ejection fraction (HFpEF). However, comorbidity profiles are usually reduced to a few comorbid disorders. Systems medicine approaches can model multiple comorbidities to help to improve our understanding of HFpEF.  

**Methods**:  
We retrospectively explored 569 comorbidities documented in 29,074 HF patients, including 8,062 HFpEF and 6,585 heart failure with reduced ejection fraction (HFrEF) patients from a German university hospital. To analyze comorbidity profiles, first we assessed differences between comorbidity profiles by HF subtype. We fit machine learning classifiers to discriminate between HFpEF and HFrEF patients to interpret important comorbidities. Then we built a network to represent comorbid relationships in HF patients (HFnet) and identified main disease clusters. Lastly, we performed novel gene candidate prediction for HFpEF by linking the HFnet to a multi-layer gene network by integrating multiple databases (HFhetnet). To corroborate HFpEF candidate genes, we collected transcriptomic data in a murine HFpEF model. We compared predicted genes with the murine disease signature as well as with the literature. 

**Results**:  
We found that comorbidity profiles vary more in association with HFpEF/HFrEF status than with age or sex. The comorbidities present in HFpEF patients were more diverse than those in HFrEF, and included neoplastic, osteologic and rheumatoid disorders, among others. Disease communities in the HFnet captured important comorbidity concepts of HF patients which could be assigned to HF subtypes, age groups and sex. Within the HFhetnet, we predicted novel gene candidates, including genes involved in fibrosis (COL3A1, MMP1, LOX), oxidative stress (NOS1, GSST1, XDH) and endoplasmic reticulum stress (ATF6). Finally, predicted genes were significantly overrepresented in the murine transcriptomic disease signature.

**Conclusions**:  
We applied systems medicine concepts to analyze comorbidity profiles in a HF patient cohort. We were able to identify disease clusters that helped to characterize HF patients. We derived a distinct comorbidity profile for HFpEF, which was leveraged to suggest novel candidate genes via network propagation. The identification of distinctive comorbidity profiles and candidate genes from routine clinical data will provide insights that may be leveraged to improve diagnosis and treatment for  HFpEF patients.

## Analysis Workflow
[Utilitary functions](https://github.com/JanLanzer/hf_comorbidity_genes/tree/master/analysis/utils)

1) Multiple Correspondence Analysis  
[MCA](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/dim_reduction/MCA_phe_reduced_script_noCMs.R)  

2) Patient Classifier  
[HFpEF_HFrEF_classifier](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/classifier/hfpef_hfref_classifier_newPID.R)  
[feature_forward_selection](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/classifier/select_estimater_cut_off.R)

3) HFnet Construction and Analysis  
[test_comorbidities_between_HFpEF/HFrEF](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/network/HFnet/perform_three_way.R)  
[construct_HFnet](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/network/HFnet/construct_HFnet_ML.R)  
[cluster_HFnet](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/network/HFnet/hfnet_clustering.R)  
[HFnet_graph_metrics](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/network/HFnet/graph_structure.R)
[Disease_network_comparison](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/network/HFnet/compare_disease_networks.R)
[Patient_cohort_and_disease_cluster_comparison](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/network/HFnet/patient_disease_cluster.R)

4) HFhetnet Construction and Analysis  
[Database_processings](https://github.com/JanLanzer/hf_comorbidity_genes/tree/master/analysis/network/data_base_processing)  
[HFhetnet_construction](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/network/HFhetnet/create_hetnet.R)  
[HFhetnet_graph_structure](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/network/HFhetnet/hetnet_characteristics.R)  
[HFhetnet_randomization](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/network/HFhetnet/randomize_layer.R)  
[LOO_CV_disease_gene_prediction](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/network/HFhetnet/validation_disease_nodes.R)
[HFhetnet_gene_prediction](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/network/HFhetnet/predict_hf_genes_ML.R)

5) RNAseq  
[processing_and_gene_annotations](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/RNAseq/process_rnaseq.R)
[Funcomics](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/RNAseq/funcomicx.R)
[Differntial_expression_analysis](https://github.com/JanLanzer/hf_comorbidity_genes/blob/master/analysis/RNAseq/DEA.R)  

## Citation
Please find our [preprint](https://doi.org/10.21203/rs.3.rs-2429581/v1)

Funded by [Informatics for Life](https://informatics4life.org/)

A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr
