import glob
import csv

#This program utilieses the output of two other programs:
#   'FullSVs_BYPatient' - A VERY non=portable program which extracts SVs from BCCA, CGI and StJude LSV TARGET datasets
#                       It's very likely that for other datasets it will be necessary to design a similar program which
#                       outputs a .csv PER PATIENT containing information on each SV. Format:
#                           line[0] - Left Chromosome 
#                           line[1] - Left Position
#                           line[2] - Right Chromosome
#                           line[3] - Right Position
#                       Any files created in this format should be portable to this program. Remember to edit paths/
#                       Naming conventions!!
#   'SVs_near_RSSs' - Much more portable, and should be usable with effectively any .bed files composed of line separated
#                     Structural Variant breakpoints.

#Folders needed: 'FullSVs_ByPatient'; 'Individual_Patient_SVs_ByRSSstatus'

#defines sources of data. Of course, if data isn't from multiple sources this list can be reduced to one entry,
#or removed completely along with 'for datasource' loops.
sources = ['CGI']

#defines default lengths used in testing for background RSS levels in the SVs_near_RSSs program. If other lengths are
#used, alter this list.
lengths = [0, 100, 300, 1000]

for datasource in sources:
    realneardict = {}
    randomneardict = {}
    
    #Creates a dictionary containing each unique, RSS-adjacent sequence in SVs_near_RSSs output.
    for nearrssfile in glob.glob(f'Unique_RSS_Site_files/{datasource}_*.csv'):
        with open(nearrssfile) as realrssfile:
            realrssreader = csv.reader(realrssfile)

            for line in realrssreader:
                chromosome = line[0]
                if chromosome == 'Chromosome':
                    pass
                else:
                    position = int(line[1]) + 45

                    realneardict[position] = chromosome

    #Does the same for each individual distance from genuine RSS in the 'lengths' list
    for distance in lengths:
        randomneardict[distance] = {}
        for randomrssfile in glob.glob(f'Random_unique_sites/{datasource}_*_{distance}tobreak_RANDOM.csv'):
            with open(randomrssfile) as simrssfile:
                randomrssreader = csv.reader(simrssfile)

                for line in randomrssreader:
                    chromosome = line[0]
                    if chromosome == 'Chromosome':
                        pass
                    else:
                        position = int(line[1]) + 45 - distance

                        randomneardict[distance][position] = chromosome

    #Creates a summary output for this program, a .csv containing information on all patients from a given datasource
    with open(f'{datasource}_SV_RSSanalysis.csv', 'w+') as fulloutfile:
        OUTfullwriter = csv.writer(fulloutfile)

        header = ['Patient', 'Diagnosis or Relapse',' ', 'Total SVs', 'SVs with both breaks near RSSs','% with both near RSSs', 'SVs with ONE break near an RSS', '% with one near an RSS','SVs with NO breaks near RSSs','% with none near RSS',' ']
        
        for distance in lengths:
            extension = [f'{distance}bp_Randoms with both breaks near RSSs',f'% {distance}bp away with both near RSSs', f'{distance}bp_Randoms with ONE break near an RSS',f'% {distance}bp away with one near an RSS', f'{distance}bp_Randoms with NO breaks near RSSs', f'% {distance}bp away with no RSSs', ' ']
            header.extend(extension)
        
        OUTfullwriter.writerow(header)
        
        #This loop cycles through each individual patient file and annotates each SV based on how many RSSs are adjacent:
        #Broken down into RSSs at one, both, or neither breakpoint
        for patientsvfile in glob.glob(f'FullSVs_ByPatient/*_{datasource}_SVs.csv'):
            caseinfo = patientsvfile.split('/')
            caseinfo = caseinfo[1].split('_')
            patientid = caseinfo[0]
            patientdorr = caseinfo[1]

            totalsvs = 0
            bothsides = 0
            oneside = 0
            neither = 0

            randombothsides = {}
            randomoneside = {}
            randomneither = {}

            for distance in lengths:
                randombothsides[distance] = 0
                randomoneside[distance] = 0
                randomneither[distance] = 0
            
            with open(patientsvfile) as fullpatientsvsfile, open(f'Individual_Patient_SVs_ByRSSstatus/{datasource}_{patientdorr}_{patientid}_SVs_near_two_RSSs.csv', 'w+') as outtwofile, open(f'Individual_Patient_SVs_ByRSSstatus/{datasource}_{patientdorr}_{patientid}_SVs_near_one_RSS.csv', 'w+') as outonefile, open(f'Individual_Patient_SVs_ByRSSstatus/{datasource}_{patientdorr}_{patientid}_SVs_near_NO_RSSs.csv', 'w+') as outneitherfile, open(f'Individual_Patient_SVs_ByRSSstatus/{datasource}_{patientdorr}_{patientid}_allSVs.csv', 'w+') as outallfile:
                individualreader = csv.reader(fullpatientsvsfile)

                #Writes an individual file for each SV in a patient with BOTH breakpoints near an RSS
                OUTtwowriter = csv.writer(outtwofile)

                #Writes an individual file for each SV in a patient with ONE breakpoint near an RSS
                OUTonewriter = csv.writer(outonefile)

                #Writes an individual file for each SV in a patient with NEITHER breakpoint near an RSS
                OUTnowriter = csv.writer(outneitherfile)

                OUTallbreaks = csv.writer(outallfile)

                for line in individualreader:
                    totalsvs += 1

                    bp1near = False
                    bp2near = False

                    patientrandom1near = {}
                    patientrandom2near = {}

                    for distance in lengths:
                        patientrandom1near[distance] = False
                        patientrandom2near[distance] = False

                    #this block of code just ensures variables brought in in a consistent format
                    bp1chr = str(line[0])
                    if "chr" not in bp1chr:
                        bp1chr = "chr" + bp1chr
                    bp1pos = int(line[1])
                    bp2chr = str(line[2])
                    if "chr" not in bp2chr:
                        bp2chr = "chr" + bp2chr
                    bp2pos = int(line[3])

                    #these try loops check whether each breakpoint can be found in the dictionary generated earlier
                    #of breaks near RSSs
                    try:
                        if realneardict[bp1pos] == bp1chr:
                            bp1near = True
                    except:
                        pass

                    try:
                        if realneardict[bp2pos] == bp2chr:
                            bp2near = True
                    except:
                        pass

                    #checks whether ONE breakpoint is near an RSS - Returns False if both or neither are near
                    oneonly = bp1near ^ bp2near

                    #Tallies the number of patient breakpoints with each proximity, and writes each SV to the appropriate OUTfile
                    if bp1near == True and bp2near == True:
                        bothsides += 1
                        OUTtwowriter.writerow(line)
                        proximity = 'Both breakpoints near RSSs'
                    elif oneonly == True:
                        oneside += 1
                        OUTonewriter.writerow(line)
                        proximity = 'One breakpoint near RSS'
                    else:
                        neither += 1
                        OUTnowriter.writerow(line)
                        proximity = 'Neither breakpoint near RSS'

                    allbreaks = line
                    allbreaks.append(proximity)
                    OUTallbreaks.writerow(allbreaks)

                    #Recreates the above analyses for each length provided, giving a background level of RSS adjacency
                    #At each distance from the genuine SV breakpoint
                    for distance in lengths:
                        try:
                            if randomneardict[distance][bp1pos] == bp1chr:
                                patientrandom1near[distance] = True
                        except:
                            pass

                        try:
                            if randomneardict[distance][bp2pos] == bp2chr:
                                patientrandom2near[distance] = True
                        except:
                            pass

                        oneonly = patientrandom1near[distance] ^ patientrandom2near[distance]

                        if patientrandom1near[distance] == True and patientrandom2near[distance] == True:
                            randombothsides[distance] += 1
                        elif oneonly == True:
                            randomoneside[distance] += 1
                        else:
                            randomneither[distance] += 1

            #Compiles a row consisting of relevant numbers for each patient and writes these to the full summary output        
            outrow = [patientid, patientdorr, ' ',totalsvs, bothsides, (bothsides/totalsvs)*100, oneside, (oneside/totalsvs)*100, neither, (neither/totalsvs)*100, ' ' ]
            
            for distance in lengths:
                extension = [randombothsides[distance], (randombothsides[distance]/totalsvs)*100, randomoneside[distance], (randomoneside[distance]/totalsvs)*100, randomneither[distance], (randomneither[distance]/totalsvs)*100, ' ']
                outrow.extend(extension)

            OUTfullwriter.writerow(outrow)    
