# A Network medicine approach to study comorbidities in heart failure with preserved ejection fraction


**Background**:
Comorbidities are expected to impact pathophysiology of heart failure (HF) with preserved ejection fraction (HFpEF). However, comorbidity profiles are usually reduced to a few comorbid disorders. Systems medicine approaches like comorbidity based gene prediction are crucial to improve our understanding of HFpEF. 
**Methods**
We retrospectively explored the comorbidities documented in 29,074 HF patients, including 8,062 HFpEF and 6,585 heart failure with reduced ejection fraction (HFrEF) patients from a German university hospital. To analyze comorbidity profiles, first we assessed differences between comorbidity profiles by HF subtype. We fit machine learning classifiers to successfully discriminate between HFpEF and HFrEF patients and interpret important comorbidities. Then we built a network to represent comorbid relationships in HF patients (HFnet) and collected unsupervised disease clusters. Lastly, we combined both approaches to perform novel gene candidate prediction for HFpEF by linkingthe HFnet to a multi-layer gene network by integrating multiple databases (HFhetnet). Within the HFhetnet we mined gene candidates by guilt by association principle. To corroborate HFpEF candidate genes we collected transcriptomic data in a murine HFpEF model and compared predicted genes with this disease signature as well as with the literature.
**REsults**: 
We found that comorbidity profiles vary to a more in association with HFpEF/HFrEF status than with age or sex. Although fewer comorbidities were recorded for HFpEF than HFrEF patients, these comorbidities were more diverse and included neoplastic, osteologic and rheumatoid disorders, among others. Disease communities in the HF comorbidity network captured important comorbidity concepts of HF patients which could be assigned to HF subtypes, age groups and sex. Within the HFhetnet, we predicted novel gene candidates, including genes involved in fibrosis (COL3A1, MMP1, LOX) and oxidative stress (NOS1, GSST1, XDH) and endoplasmic reticulum stress (ATF6). Finally, we corroborated some predicted genes in myocardial gene expression data. 
In summary, we applied systems medicine concepts to analyze comorbidity profiles in a multimorbid HF patient cohort. We were able to identify the main disease groups of HF comorbidities, and distinct comorbidity profiles for HFpEF which were leveraged to suggest novel candidate genes via network propagation. The identification of distinctive comorbidity profiles as well as candidate genes from routine clinical data will provide insights that may be leveraged to improve diagnosis and treatment for  HFpEF patients.




**Analysis Workflow**

Processing, QC and integration
1) sample wise QC [run_samplewise_processing.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/sample_integration/run_sample_wise_preprocessing.R) 
2) first sample integration [harmony_integration_across_sample.R](https://github.com/saezlab/scell_hfpef/blob/main/analysis/sample_integration/harmony_integration_across_sample.R)



A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr
