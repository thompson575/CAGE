## CAGE: Classification After Gene Expression

### Project Aim

CAGE explores ways of using data from a gene expression study for subject classification. In particular, it investigates the pros and cons of basing the classification on the principal components of the gene expression.

### Example

The gene expression study used in the project comes from,

Xu K, Shi X, Husted C, Hong R et al.  
**Smoking modulates different secretory subpopulations expressing SARS-CoV-2 entry genes in the nasal and bronchial airways.**  
Sci Rep 2022 Oct 28;12(1):18168. 

The authors made their data available on the GEO archive as GSE210271. The CAGE analysis is based on the Series Matrix File downloaded from <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210271>.  

On GEO the study design is described as,

"mRNA gene expression from 505 nasal epithelial brushings was profiled using Affymetrix Gene 1.0 ST microarrays. Samples were collected from patients in the Airway Epithelial Gene Expression in the Diagnosis of Lung Cancer (AEGIS) trials (AEGIS-1 and AEGIS-2), two independent, prospective, multicenter, observational studies. 375 nasal samples (243 with lung cancer, 132 with benign lung disease) were patients in the AEGIS-1 trial and 130 nasal samples (66 with lung cancer, 64 with benign lung disease) were from patients in the AEGIS-2 trial."  

The classification problem is to use the gene expression measurements to distinguish lung cancer from benign lung disease.

