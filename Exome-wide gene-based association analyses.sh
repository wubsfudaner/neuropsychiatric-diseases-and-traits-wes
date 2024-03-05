#!/bin/bash
pheno="pheno"
conda activate SAIGE

#for binary phenotype

for i in {1..22};
do
Rscript step1_fitNULLGLMM.R \
    --sparseGRMFile=UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --plinkFile=ukb_wes_chr${i}_sample_qc_final_unrelated \
    --useSparseGRMtoFitNULL=TRUE \
    --phenoFile=${pheno}_wes_data.csv \
    --phenoCol=${pheno} \
    --covarColList=age_onset,sex,batch,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --qCovarColList=sex,batch  \
    --sampleIDColinphenoFile=eid \
    --traitType=binary \
    --isCateVarianceRatio=TRUE \
    --nThreads=20 \
    --LOCO=TRUE \
    --SampleIDIncludeFile=british_id_sample.txt \
    --outputPrefix=/STEP1/${pheno}_chr${i} \
    --IsOverwriteVarianceRatioFile=TRUE

#Exome-wide gene-based association analyses    
Rscript step2_SPAtests.R \
    --bedFile=ukb_wes_chr${i}_sample_qc_final_unrelated.bed \
    --bimFile=ukb_wes_chr${i}_sample_qc_final_unrelated.bim \
    --famFile=ukb_wes_chr${i}_sample_qc_final_unrelated.fam \
    --SAIGEOutputFile=/STEP2/${pheno}_chr${i} \
    --AlleleOrder=alt-first \
    --minMAF=0 \
    --minMAC=0.5 \
    --GMMATmodelFile=/STEP1/${pheno}_chr${i}.rda \
    --varianceRatioFile=/STEP1/${pheno}_chr${i}.varianceRatio.txt \
    --sparseGRMFile=UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --groupFile=SnpEff_gene_group_chr${i}.txt \
    --annotation_in_groupTest="lof,missense,missense:lof" \
    --maxMAF_in_groupTest=0.000005,0.00001,0.0001,0.001,0.01 \
    --is_output_markerList_in_groupTest=TRUE \
    --LOCO=F \
    --is_fastTest=TRUE
done


#for quantitative traits

for i in {1..22};
do
Rscript step1_fitNULLGLMM.R \
    --sparseGRMFile=UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --plinkFile=ukb_wes_chr${i}_sample_qc_final_unrelated \
    --useSparseGRMtoFitNULL=TRUE \
    --phenoFile=${pheno}_wes_data.csv \
    --phenoCol=${pheno} \
    --covarColList=age,sex,batch,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --qCovarColList=sex,batch  \
    --sampleIDColinphenoFile=eid \
    --traitType=quantitative \
    --invNormalize=T \
    --isCateVarianceRatio=F \
    --nThreads=20 \
    --LOCO=TRUE \
    --SampleIDIncludeFile=british_id_sample.txt \
    --outputPrefix=/STEP1/${pheno}_chr${i} \
    --IsOverwriteVarianceRatioFile=TRUE

#Exome-wide gene-based association analyses    
Rscript step2_SPAtests.R \
    --bedFile=ukb_wes_chr${i}_sample_qc_final_unrelated.bed \
    --bimFile=ukb_wes_chr${i}_sample_qc_final_unrelated.bim \
    --famFile=ukb_wes_chr${i}_sample_qc_final_unrelated.fam \
    --SAIGEOutputFile=/STEP2/${pheno}_chr${i} \
    --AlleleOrder=alt-first \
    --minMAF=0 \
    --minMAC=0.5 \
    --GMMATmodelFile=/STEP1/${pheno}_chr${i}.rda \
    --varianceRatioFile=/STEP1/${pheno}_chr${i}.varianceRatio.txt \
    --sparseGRMFile=UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --groupFile=SnpEff_gene_group_chr${i}.txt \
    --annotation_in_groupTest="lof,missense,missense:lof" \
    --maxMAF_in_groupTest=0.000005,0.00001,0.0001,0.001,0.01 \
    --is_output_markerList_in_groupTest=TRUE \
    --LOCO=F \
    --is_fastTest=TRUE
done