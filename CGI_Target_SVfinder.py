import statistics
import csv

#This program was designed VERY SPECIFICALLY to pull out SV breakpoints from the input .csv file
#It is unlikely to be broadly applicable

#Inputs: CGI_SV_Breakpoints.ssv file, from TARGET
#       Target_Data_Matrix.csv from TARGET
#       DESeq_Metadata_Expanded file (from matrix_for_deseq2/metadata_expander)

#File outputs: a Metadata file containing information on the number of SVs by various categories
#               (ie subtype, relapse vs diagnosis) (.txt)
#               A full list of breakpoints belonging to annotated cases (.csv)
#               A full list of breakpoints belonging to unannotated cases (.csv)

def main():
    with open('CGI_Target/CGI_SV_Breakpoints.csv') as inputfile, open ('DESeq_Metadata_expanded.csv') as metafile, open('CGI_Target/CGI_Target_SVmetadata.txt', 'w+') as svmetadata, open ('CGI_Target/CGI_Target_SVBreakpoints.csv', 'w+') as breakpoints, open('CGI_Target/CGI_Target_UNCLASSIFIED_breakpoints.csv', 'w+') as unclasssified:
        
        #First we prep our I/O files

        inreader = csv.reader(inputfile)
        next(inreader)

        metareader = csv.reader(metafile)
        next(metareader)

        OUTbreaks = csv.writer(breakpoints)
        OUTunclassified = csv.writer(unclasssified)

        #OUT files header written
        breakshead = ['Patient ID', 'Subtype', 'Diagnosis or Relapse', 'Chromosome', 'Position', 'Nearby Gene']
        OUTbreaks.writerow(breakshead)

        unclassifiedhead = ['Patient ID', 'Chromosome', 'Position', 'Nearby Gene']
        OUTunclassified.writerow(unclassifiedhead)

        #Various necessary variables defined
        patients = []
        nobreakslist = []
        shortidlist =[]
        annotatedpatients = []
        unannotatedpatients = []
        subtypebypatient = {}
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

        #Each case annotated with subtype based on DESeq Metadata
        for line in metareader:
            patientid = line[0]
            patientidshort = patientid[:-8]

            subtypebypatient[patientidshort] = line[1]

        for line in inreader:
            totalbreakpoints += 2

            #Variables found in .csv file
            patientid = line[0]
            patientidshort = str(patientid[:-8])
            chrleft = line[1]
            posleft = line[2]
            chrright = line[5]
            posright = line[6]
            geneleft = line[18]
            generight = line[19]

            #Number of total breaks per patient incremented fo each found SV
            if patientid not in patients:
                breaksbypatient[patientid] = 2
                patients.append(patientid)
            else:
                breaksbypatient[patientid] += 2
            
            if patientidshort not in shortidlist:
                shortidlist.append(patientidshort)
            
            #Checks RNAseq sample numbers against data matrix to get diagnosis/relapse info
            #SVs from diagnosis/relapse are counted for Metadata file
            with open('TARGET_Data_Matrix.csv') as matrix:
                matrixreader = csv.reader(matrix)

                for line in matrixreader:
                    targetpatient = str(line[0])

                    if targetpatient == patientidshort:

                        if patientid in line[31]:
                            dorrbyid = 'diagnosis'
                            possibleids = line[5].split(',')

                            if patientidshort not in diagnosis:
                                diagnosis.append(patientidshort)
                                breaksbypatientdiag[patientidshort] = 2

                            else:
                                breaksbypatientdiag[patientidshort] += 2

                            break

                        elif patientid in line[33]:
                            dorrbyid = 'relapse'
                            possibleids = line[6].split(',')

                            if patientidshort not in relapse:
                                relapse.append(patientidshort)
                                breaksbypatientrelapse[patientidshort] = 2

                            else:
                                breaksbypatientrelapse[patientidshort] += 2

                            break

                    else:
                        dorrbyid = 'N/A - not in sample matrix'    

            #This try/except checks for known subtype info for each patient
            try:
                #All annotated breakpoints are counted        
                subtype = subtypebypatient[patientidshort]
                annotatedbreaks += 2
                outlineleft = [patientid, subtype, dorrbyid, chrleft, posleft, geneleft]
                outlineright = [patientid, subtype, dorrbyid, chrright, posright, generight]

                #Writes breakpoints to the annotated SV breakpoint file
                OUTbreaks.writerow(outlineleft)
                OUTbreaks.writerow(outlineright)

                if patientid not in annotatedpatients:
                    annotatedpatients.append(patientid)

            #Unannotated breaks are all handled seperately, in a similar way
            except: 
                unannotatedbreaks += 2

                if patientid not in unannotatedpatients:
                    unannotatedpatients.append(patientid)

                unclassoutleft = [patientid, chrleft, posleft, geneleft]
                unclassoutright = [patientid, chrright, posright, generight]

                OUTunclassified.writerow(unclassoutleft)
                OUTunclassified.writerow(unclassoutright)

        #Generates a list of breaks present in each patient, at diagnosis vs relapse
        for patient in patients:
            nobreakslist.append(breaksbypatient[patient])

        for patient in shortidlist:
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
        svmetadata.write(f'{len(patients)} samples \n')
        svmetadata.write(f'Samples from {len(shortidlist)} patients \n')
        svmetadata.write(f'{len(annotatedpatients)} samples had known subtypes \n')
        svmetadata.write(f'{len(unannotatedpatients)} samples were unclassified - NO AVAILABLE RNAseq data \n')

        svmetadata.write('\n ANNOTATED CASES \n')
        for patient in annotatedpatients:
            svmetadata.write(f'{patient} \n')

        svmetadata.write('\n UNANNOTATED CASES \n')
        for patient in unannotatedpatients:
            svmetadata.write(f'{patient} \n')

        svmetadata.write(f'\n {len(diagnosis)} samples from diagnosis \n')
        svmetadata.write(f'{len(relapse)} samples from relapse\n')

        svmetadata.write('\n Diagnosis Samples \n')
        for patient in diagnosis:
            svmetadata.write(f'{patient} \n')

        svmetadata.write('\n Relapse Samples \n')
        for patient in relapse:
            svmetadata.write(f'{patient} \n')

        svmetadata.write(f'\n {totalbreakpoints} total breakpoints \n')
        svmetadata.write(f'{annotatedbreaks} breakpoints from annotated cases \n')
        svmetadata.write(f'{unannotatedbreaks} breakpoints from unannotated cases \n')
        svmetadata.write(f'{avgbreaks} breakpoints on average per patient \n')
        svmetadata.write(f'{medbreaks} median breaks \n')
        svmetadata.write(f'{avgdiagbreaks} breakpoints on average per DIAGNOSIS sample\n')
        svmetadata.write(f'{meddiagbreaks} median breaks in DIAGNOSIS samples \n')
        svmetadata.write(f'{avgrelapsebreaks} breakpoints on average per RELAPSE sample \n')
        svmetadata.write(f'{medrelapsebreaks} median breaks in RELAPSE samples \n \n')
        svmetadata.write('no. breakpoints by patient \n \n')

        for patient in patients:
            svmetadata.write(f'{patient}:      {breaksbypatient[patient]} \n')
        
        svmetadata.write('\n \n no. breakpoints by patient, split diagnosis v relapse \n \n ')

        for patient in shortidlist:
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

        
