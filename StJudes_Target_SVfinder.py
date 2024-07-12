import csv
import statistics
import sys

#This program was designed VERY SPECIFICALLY to pull out SV breakpoints from the input .tsv file
#It is unlikely to be broadly applicable

#Inputs: St Judes WGS_SV.tsv file, from TARGET
#       DESeq_Metadata_Expanded file (from matrix_for_deseq2/metadata_expander)

#File outputs: a Metadata file containing information on the number of SVs by various categories
#               (ie subtype, relapse vs diagnosis) (.txt)
#               A full list of breakpoints belonging to annotated cases (.csv)
#               A full list of breakpoints belonging to unannotated cases (.csv)

def main():
    with open(sys.argv[1]) as inputfile, open (sys.argv[2]) as metafile, open('StJudes_Target/StJudes_Target_SVmetadata.txt', 'w+') as svmetadata, open ('StJudes_Target/StJudes_Target_SVBreakpoints.csv', 'w+') as breakpoints, open('StJudes_Target/StJudes_Target_UNCLASSIFIED_breakpoints.csv', 'w+') as unclasssified:
        
        #First we prep our I/O files
        inreader = csv.reader(inputfile, delimiter = '\t')
        next(inreader)

        metareader = csv.reader(metafile)
        next(metareader)

        OUTbreaks = csv.writer(breakpoints)
        OUTunclassified = csv.writer(unclasssified)

        #OUT files header written
        breakshead = ['Patient ID', 'Subtype', 'Diagnosis or Relapse', 'Chromosome', 'Position', 'Nearby Gene']
        OUTbreaks.writerow(breakshead)

        unclassifiedhead = ['Patient ID', 'Diagnosis or Relapse', 'Chromosome', 'Position', 'Nearby Gene']
        OUTunclassified.writerow(unclassifiedhead)
        
        #Various necessary variables defined
        patients = []
        nobreakslist = []
        annotatedpatients = []
        unannotatedpatients = []
        subtypebypatient = {}
        DorRbypatient = {}
        breaksbypatient = {}
        diagnosis = []
        relapse = []
        totalbreakpoints = 0
        annotatedbreaks = 0
        unannotatedbreaks = 0
        diagnobreakslist = []
        relapsenobreakslist = []
        breaksbypatientdiag = {}
        breaksbypatientrelapse = {}

        #Each case annotated with subtype and diagnosis/relapse based on DESeq Metadata
        for line in metareader:
            patientid = line[0]
            patientid = patientid[:-8]
            patientid = patientid[7:]
            subtypebypatient[patientid] = line[1]
            DorRbypatient[patientid] = line[2]

        for line in inreader:
            totalbreakpoints += 2

            #Variables found in .tsv file
            patientid = line[1]
            dorr = line[2]
            chrleft = line[5]
            posleft = line[6]
            chrright = line[10]
            posright = line[11]
            geneleft = line[4]
            generight = line[9]

            #Number of total breaks per patient incremented fo each found SV
            if patientid not in patients:
                patients.append(patientid)
                breaksbypatient[patientid] = 2

            else:
                breaksbypatient[patientid] += 2

            #This try/except checks for known subtype info for each patient
            try:
                #All annotated breakpoints are counted
                subtype = subtypebypatient[patientid]
                annotatedbreaks += 2

                if patientid not in annotatedpatients:
                    annotatedpatients.append(patientid)

                #checks and increments breaks from diagnosis
                if dorr == 'diagnosis':

                    if patientid not in diagnosis:
                        diagnosis.append(patientid)
                        breaksbypatientdiag[patientid] = 2

                    else:
                        breaksbypatientdiag[patientid] += 2

                #checks and increments breaks from relapse
                elif dorr == 'relapse':
                    if patientid not in relapse:
                        relapse.append(patientid)
                        breaksbypatientrelapse[patientid] = 2

                    else:
                        breaksbypatientrelapse[patientid] += 2

                #Writes breakpoints to the annotated SV breakpoint file
                outlineleft = [patientid, subtype, dorr, chrleft, posleft, geneleft]
                outlineright = [patientid, subtype, dorr, chrright, posright, generight]
                OUTbreaks.writerow(outlineleft)
                OUTbreaks.writerow(outlineright)

            #Unannotated breaks are all handled seperately, in a similar way
            except: 
                unannotatedbreaks += 2

                if patientid not in unannotatedpatients:
                    unannotatedpatients.append(patientid)

                if dorr == 'diagnosis':

                    if patientid not in diagnosis:
                        diagnosis.append(patientid)
                        breaksbypatientdiag[patientid] = 2

                    else:
                        breaksbypatientdiag[patientid] += 2


                elif dorr == 'relapse':

                    if patientid not in relapse:
                        relapse.append(patientid)
                        breaksbypatientrelapse[patientid] = 2

                    else:
                        breaksbypatientrelapse[patientid] += 2

                unclassoutleft = [patientid, dorr, chrleft, posleft, geneleft]
                unclassoutright = [patientid, dorr, chrright, posright, generight]

                OUTunclassified.writerow(unclassoutleft)
                OUTunclassified.writerow(unclassoutright)

        #Generates a list of breaks present in each patient at diagnosis vs relapse
        for patient in patients:
            nobreakslist.append(breaksbypatient[patient])

            try:
                diagnobreakslist.append(breaksbypatientdiag[patient])
            except:
                pass

            try:
                relapsenobreakslist.append(breaksbypatientrelapse[patient])
            except:
                pass

        #calculates some statsfor the metadata file
        avgbreaks = (sum(nobreakslist)/len(nobreakslist))
        medbreaks = statistics.median(nobreakslist)

        avgdiagbreaks = (sum(diagnobreakslist)/len(diagnobreakslist))
        meddiagbreaks = statistics.median(diagnobreakslist)

        avgrelapsebreaks = (sum(relapsenobreakslist)/len(relapsenobreakslist))
        medrelapsebreaks = statistics.median(relapsenobreakslist)

        #Everything from here on down writes pertinent info to the Metadata file
        svmetadata.write(f'{len(patients)} patients in list \n')
        svmetadata.write(f'{len(annotatedpatients)} had known subtypes \n')
        svmetadata.write(f'{len(unannotatedpatients)} were unclassified - NO AVAILABLE RNAseq data \n')
        svmetadata.write(f'{len(diagnosis)} samples from diagnosis \n')
        svmetadata.write(f'{len(relapse)} samples from relapse\n')
        svmetadata.write(f'{totalbreakpoints} total breakpoints \n')
        svmetadata.write(f'{annotatedbreaks} breakpoints from annotated cases \n')
        svmetadata.write(f'{unannotatedbreaks} breakpoints from unannotated cases \n')

        svmetadata.write(f'{avgbreaks} breakpoints on average per patient \n')
        svmetadata.write(f'{medbreaks} median breaks \n')
        svmetadata.write(f'{avgdiagbreaks} breakpoints on average per DIAGNOSIS sample\n')
        svmetadata.write(f'{meddiagbreaks} median breaks in DAGNOSIS samples \n')
        svmetadata.write(f'{avgrelapsebreaks} breakpoints on average per RELAPSE sample \n')
        svmetadata.write(f'{medrelapsebreaks} median breaks in RELAPSE samples \n \n')
        svmetadata.write('no. breakpoints by patient \n \n')

        for patient in patients:
            svmetadata.write(f'{patient}:      {breaksbypatient[patient]} \n')
        
        svmetadata.write('\n \n no. breakpoints by patient, split diagnosis v relapse \n \n ')

        for patient in patients:
            try:
                diagnosisbreaks = breaksbypatientdiag[patient]
            except:
                diagnosisbreaks = 'N/A'
            
            try:
                relapsebreaks = breaksbypatientrelapse[patient] 
            except:
                relapsebreaks = 'N/A'
            
            svmetadata.write(f'{patient}:       {diagnosisbreaks} breaks at diagnosis       {relapsebreaks} breaks at relapse \n')

if __name__=='__main__':
    main()