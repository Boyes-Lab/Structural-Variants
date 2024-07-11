import glob
import csv

#This program is VERY specifically designed for the specifics of these datasets. Broad applicability is unlikely.

#This program outputs one .csv file per patient, one FULL SV (ie two breakpoints) per line
#These are output to a folder in the current working directory, "FullSVs_ByPatient", which must be created to run the program

#This dictionary stored ALL patient SV information
fullpatientdict = {}

#SVs from the StJudes SV.tsv file are collected
#This file INCLUDES whether a case is diagnosis or relapse
with open('StJudes_Target/TARGET_BALL_WGS_SV.tsv') as stjudesfile:
    stjudereader = csv.reader(stjudesfile, delimiter = '\t')
    next(stjudereader)

    for line in stjudereader:
        patientid = line[1]
        dorr = line[2].upper()
        updatedid = patientid + '_' + dorr + '_StJude'

        if updatedid not in fullpatientdict:
            fullpatientdict[updatedid] = []

        bp1chr = line[5]
        bp1pos = line[6]
        bp2chr = line[10]
        bp2pos = line[11]

        fullsv = [bp1chr, bp1pos, bp2chr, bp2pos]

        fullpatientdict[updatedid].append(fullsv)

#The CGI annotated breaks file (from CGI_Target_SVFinder) contains information on whether
#  a given sample is diagnosis or relapse
CGIannotations = {}

with open('OUTPUTS/CGI_annotated_breaks.csv') as cgiannofile:
    cgiannoreader = csv.reader(cgiannofile)

    for line in cgiannoreader:
        patientid = line[0]
        dorr = line[2]
        if 'N/A' in dorr:
            dorr = 'UNANNOTATED'
        if patientid not in CGIannotations:
            CGIannotations[patientid] = dorr

#SVs from the CGI breakpoints.csv file are collected and annotated as diagnosis, relapse or unannotated
with open('CGI_Target/CGI_SV_Breakpoints.csv') as cgifile:
    cgireader = csv.reader(cgifile)

    for line in cgireader:
        patientid = line[0]

        try:
            dorr = CGIannotations[patientid]
        except:
            dorr = 'UNANNOTATED'
        
        updatedid = patientid + '_' + dorr + '_CGI'

        if updatedid not in fullpatientdict:
            fullpatientdict[updatedid] = []

        bp1chr = line[1]
        bp1pos = line[2]
        bp2chr = line[5]
        bp2pos = line[6]

        fullsv = [bp1chr, bp1pos, bp2chr, bp2pos]

        fullpatientdict[updatedid].append(fullsv)

#SVs from individual patient BCCA somatic breakpoints.tsv file are collected
#Note that ALL BCCA cases are from diagnosis only
for patienttsv in sorted(glob.glob(f'BCCA_Target/BCCA_Large_Somatic/*.tsv')):

    patientid = patienttsv.split('/')
    patientid = patientid[2].split('.')
    patientid = patientid[0]
    updatedid = patientid +'_DIAGNOSIS_BCCA'

    if updatedid not in fullpatientdict:
        fullpatientdict[updatedid] = []

    with open(patienttsv) as inputfile:
        inreader = csv.reader(inputfile, delimiter = '\t')
        next(inreader)

        for line in inreader:
            workingline = line[0].split('_')
            geneleft = workingline[2]
            generight = workingline[3]
            leftright = workingline[1].split('|')
            left = leftright[0].split(':')
            bp1chr = left[0]
            bp1pos = left[1]
            right = leftright[1].split(':')
            bp2chr = right[0]
            bp2pos = right[1]

            fullsv = [bp1chr, bp1pos, bp2chr, bp2pos]

            fullpatientdict[updatedid].append(fullsv)

#Finally, an output file is written for all SVs from a single patient
for patient in fullpatientdict:
    with open(f'FullSVs_ByPatient/{patient}_SVs.csv', 'w+') as output:
        outwriter = csv.writer(output)

        for structuralvariant in fullpatientdict[patient]:  
            outwriter.writerow(structuralvariant)