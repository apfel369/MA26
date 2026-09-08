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
Finaler Ablauf:
- HBeeID_Multi-Sample-VCF_based-on-lifted-excel-snps_AmelHAv31.sh
- v2-filter-multi-sample-vcf-amelhav31.py
- liftover_to_amel45-and-match_ID.sh
- v2-indels-entfernen-hbeeid-input-data.sh
- v2-genotypmatrix-hbeeid-erzeugen.py

Durchlauf F:
- HBeeID-272-amel45-snps_multi-sample-VCF.sh
- do-biallelic_indels-to-NA-genotype_amel45.sh
- v2-filter-multi-sample-vcf-amelhav31.py wurde analog genutzt (mit angepassten Pfaden und angepassten Dateien) 
- v2-genotypmatrix-hbeeid-erzeugen-amel45.py

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

 ### 4.4.3 Admixture-Analyse ###
 - NGSadmix_code_MAF-5.txt
