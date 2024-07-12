import statistics
import csv
import glob

#This program was designed VERY SPECIFICALLY to pull out SV breakpoints from input .tsv files - BCCA separates these
#BY PATIENT, and also into large vs short SVs

#Inputs: somatic.large.summary.tsv file PER PATIENT from TARGET, in BCCA_TARGET/BCCA_Large_Somatic folder
#        somatic.short.summary.tsv file PER PATIENT from TARGET, in BCCA_TARGET/BCCA_Short_Somatic folder
#       DESeq_Metadata_Expanded file (from matrix_for_deseq2/metadata_expander)

#File outputs: a Metadata file containing information on the number of SVs by various categories
#               (ie subtype, relapse vs diagnosis) (.txt)
#               A full list of breakpoints belonging to annotated cases (.csv)
#               A full list of breakpoints belonging to unannotated cases (.csv)
#                   - Two of each of these files (one for Large, one for Short)

def main():
    lengths = ('Large', 'Short')

    for length in lengths:
        with open(f'BCCA_Target/BCCA_Target_{length}_SVmetadata.txt', 'w+') as svmetadata, open (f'BCCA_Target/BCCA_Target_{length}_SVBreakpoints.csv', 'w+') as breakpoints, open(f'BCCA_Target/BCCA_Target_{length}_UNCLASSIFIED_breakpoints.csv', 'w+') as unclasssified:
            #First we prep our OUT files
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

            #In this case, input files are opened later to combine all breaks from seperate files into one.
            for patienttsv in sorted(glob.glob(f'BCCA_Target/BCCA_{length}_Somatic/*.tsv')):

                with open(patienttsv) as inputfile, open ('DESeq_Metadata_expanded.csv') as metafile: 
                    #Now we prep the input files
                    inreader = csv.reader(inputfile, delimiter = '\t')
                    next(inreader)

                    metareader = csv.reader(metafile)
                    next(metareader)

                    #Each case annotated with subtype and diagnosis/relapse based on DESeq Metadata
                    for line in metareader:
                        patientid = line[0]
                        patientidshort = patientid[:-8]

                        subtypebypatient[patientidshort] = line[1]
                        DorRbypatient[patientidshort] = line[2]

                    #gets patient ID from file name
                    currenttsvid = patienttsv.split('/')
                    currenttsvid = currenttsvid[2]
                    currenttsvid = currenttsvid[:16]

                    for line in inreader:

                        #variables found in .tsv file
                        workingline = line[0].split('_')
                        geneleft = workingline[2]
                        generight = workingline[3]
                        leftright = workingline[1].split('|')
                        left = leftright[0].split(':')
                        chrleft = left[0]
                        posleft = left[1]
                        right = leftright[1].split(':')
                        chrright = right[0]
                        posright = right[1]

                        #Number of total breaks per patient incremented fo each found SV
                        if currenttsvid not in patients:
                            breaksbypatient[currenttsvid] = 2
                            patients.append(currenttsvid)
                        else:
                            breaksbypatient[currenttsvid] += 2

                        #This try/except checks for known subtype info for each patient
                        try:      
                            subtype = subtypebypatient[currenttsvid]
                            dorrbyid = DorRbypatient[currenttsvid]

                            #checks and increments breaks from diagnosis
                            if dorrbyid == 'Diagnosis':

                                if currenttsvid not in diagnosis:
                                    diagnosis.append(currenttsvid)

                                breaksbypatientdiag[currenttsvid] = breaksbypatient[currenttsvid]

                            #checks and increments breaks from relapse
                            elif dorrbyid == 'Relapse':

                                if currenttsvid not in relapse:
                                    relapse.append(currenttsvid)

                                breaksbypatientrelapse[currenttsvid] = breaksbypatient[currenttsvid]

                            #All annotated breakpoints are counted  
                            annotatedbreaks += 2

                            outlineleft = [currenttsvid, subtype, dorrbyid, chrleft, posleft, geneleft]
                            outlineright = [currenttsvid, subtype, dorrbyid, chrright, posright, generight]

                            OUTbreaks.writerow(outlineleft)
                            OUTbreaks.writerow(outlineright)

                            #Writes breakpoints to the annotated SV breakpoint file
                            if currenttsvid not in annotatedpatients:
                                annotatedpatients.append(currenttsvid)

                        #Unannotated breaks are all handled seperately, in a similar way
                        except: 
                            unannotatedbreaks += 2

                            if patientid not in unannotatedpatients:
                                unannotatedpatients.append(patientid)

                            unclassoutleft = [currenttsvid, chrleft, posleft, geneleft]
                            unclassoutright = [currenttsvid, chrright, posright, generight]

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

            totalbreakpoints = sum(nobreakslist)

            #calculates some statsfor the metadata file
            avgbreaks = (sum(nobreakslist)/len(nobreakslist))
            medbreaks = statistics.median(nobreakslist)

            try:
                avgdiagbreaks = (sum(diagnobreakslist)/len(diagnobreakslist))
                meddiagbreaks = statistics.median(diagnobreakslist)
            except:
                avgdiagbreaks = 'Cannot be calculated'
                meddiagbreaks = 'Cannot be calculated'

            try:
                avgrelapsebreaks = (sum(relapsenobreakslist)/len(relapsenobreakslist))
                medrelapsebreaks = statistics.median(relapsenobreakslist)
            except:
                avgrelapsebreaks = 'Cannot be calculated'
                medrelapsebreaks = 'Cannot be calculated'

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

            svmetadata.write(f'\n{totalbreakpoints} total breakpoints \n')
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