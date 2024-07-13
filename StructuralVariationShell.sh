#!/bin/bash

#Inputs as follows:
#Datasource - one of BCCA, StJudes or CGI
datasource=$1
#The filename (or prefix) of your data downloaded from TARGET. 
    #For BCCA, "Somatic.large.summary" and "Somatic.small.summary", with identical prefixes, FOR EACH PATIENT
    #For CGI, 
    #For St Judes, the "TARGET_BALL_WGS_SV" file
dataset_of_interest=$2
#A metadata file describing the datasets in question 
metadata=$3
#The TARGET data matrix, defining patient characteristics (not provided due to patient information security)
target_matrix=$4
#A folder containing all patient.bed files 
bedfolder=$5
#The genome file being used (in FASTA format)
genome=$6

#First we find Structural Variants in a datasource-specific manner

mkdir OUTPUTS
mkdir Individual_Patient_SVs_ByRSSstatus

if [[ "$datasource" == "CGI" ]]; then
    mkdir CGI_Target
    python CGI_target_SVfinder.py $dataset_of_interest $metadata $target_matrix
    svfile="OUTPUTS/CGI_SV_Breakpoints"
elif [[ "$datasource" == "StJude" ]]; then
    mkdir BCCA_Target
    python BCCA_target_SVfinder.py $dataset_of_interest $metadata 
    svfile="NA"
elif [[ "$datasource" == "BCCA" ]]; then
    mkdir StJudes_Target
    python CGI_target_SVfinder.py $dataset_of_interest $metadata $target_matrix
    svfile="NA"
fi

#Then we get the full complement of Structural Variants for each patient
mkdir FullSVs_ByPatient
python FullSVs_ByPatient.py $datasource $dataset_of_interest $svfile

#Third, we find which of these are near RSS sites using RSSsite
python SVs_near_RSSs.py -g $genome -b $bedfolder

#Finally, we produce some outputs summarising SV content and proximity to RSS by patient
python Finalise_SVs.py $datasource



