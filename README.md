Nachfolgend ist die Zuordnung der genutzten Skripte zu den jeweilig erfolgten Methoden. 

## 4.1 Daten-Vorverarbeitung ##
- fastp.sh
- mapping_bowtie2_buckfast-and-upset-samples.sh
- Picard_mark-dupes.sh

## 4.2 Varianten-Identifikation ##
- variant-calling_single-sample.sh
- variant-filtering_single-sample-VCFs.sh

## 4.3 Quantitativer SNV-Vergleich ##
Die dazu gekommenen SRA-Daten wurden zunächst wie folgt verarbeitet: 
- fastp.sh
- mapping_bowtie2_buckfast-and-upset-samples.sh
- Picard_mark-dupes.sh
- variant-calling_single-sample.sh
- variant-filtering_single-sample-VCFs.sh

Mit Buckfast dann gemeinsam:
- Remove-InDels.txt
- ComplexUpset.R

## 4.4 Populationsstrukturanalyse ##
### 4.4.1 Klassifikation ###
Skripte sind unter anderen Branch gelistet

 ### 4.4.2 PCA ###
Die neuen Referenzproben Monticola und vom CNP-Paper zunächst wie folgt verarbeitet: 
- fastp.sh
- mapping_bowtie2_buckfast-and-upset-samples.sh
- Picard_mark-dupes.sh

#### 4.4.2.1 Genotypbasierter Ansatz ####
- variant-calling_multi-sample-VCF.sh
- Remove-InDels.txt
- vcf-filter_GT-workflow_step-01.sh
- vcf-filter_GT-workflow_step-02.txt
- LD-and-PCA_plink2.txt

#### 4.4.2.2 Likelihood-basierter Ansatz ####
- GL-Pipeline.sh
