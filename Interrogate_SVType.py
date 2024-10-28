import csv
import glob
import statistics

from collections import defaultdict

#Folders needed: 'FullSVs_DNADamageExpanded'; 'FullSVs_DNADamageExpanded/SplitByRSSProximity'; 
#                'FullSVs_DNADamageExpanded/SplitByRSSProximity/ByPatientComplied'
#                'FullSVs_DNADamageExpanded/SplitByRSSProximity/Summaries'


#This function checks whether a given genomic sequence corresponds to an Ig/Tcr locus
#Co ordinates based on hg19 alignment - CHANGE THESE AS APPROPRIATE FOR ALTERNATE GENOME BUILDS
#Notably, this is very obviously easily modified to search for ANY locus of interest
def locifind(position, chromosome):

    #Regions are listed in format:
        # (Region name, chromosome, chromosomal region)
    regions = [
        ("IgH", "chr14", range(106052274,107388051) ),
        ('TrA or TrD', "chr14", range(22090057,23021075)),
        ("TrB", 'chr7', range(141998851,142510972)),
        ("TrC", 'chr7', range(38279625,38407656)),
        ("TrC", 'chr7', range(38279625,38407656)),
        ("IgK", 'chr2', range(89156874,90274235) ),
        ("IgL", 'chr22',range(22380474,23265085))
    ]

    for locus, chrom_target, region_target in regions:
        if (chromosome == chrom_target) and (position in region_target):
            return locus

    return 'No match to Ig or Tcr regions'

#VERY simple function to work out value1 as a percentage of value2
def percentage(value1, value2):
    return((value1/value2)*100)

#Taking in the 'Unique_RSS_Site_files" folder from SVs_near_RSSs.py, gathers a dictionary with the format:
#   datasource[chromosomal position] = chromosome
def gather_RSS_proximity(sources, proximity_folder):
    realneardict = defaultdict(dict)

    for datasource in sources:

        for nearrssfile in glob.glob(f'{proximity_folder}/{datasource}_*.csv'):
            
            with open(nearrssfile) as realrssfile:
                realrssreader = csv.reader(realrssfile)
                
                for line in realrssreader:
                    chromosome = line[0]

                    #This just skips over the Header line
                    if chromosome == 'Chromosome':
                        pass
                    
                    #This corrects for the window applied in the SVs_near_RSSs program, and creates a dictionary
                    #of all breakpoints which are near RSSs
                    else:
                        position = int(line[1]) + 45

                        realneardict[datasource][position] = chromosome
    
    return realneardict

def main():

    #defines sources of data. Of course, if data isn't from multiple sources this list can be reduced to one entry,
    #or removed completely along with 'for datasource' loops.
    sources = ['StJude', 'BCCA', 'CGI']

    #defines SVs by their proximity to RSS sites, as identified using 'SVs_near_RSSs_complete.py' program. 
    #This program is assumed to be the source of all input data
    proximity= ['near_NO_RSSs', 'near_one_RSS', 'near_two_RSSs']

    #Needed for the production of metadata files on proportions of inter- vs intra-chromosomal events
    chromosomalchar = ['Interchromosomal', 'Intrachromosomal']


    individualpatientdict = {}

    locicountdict = defaultdict(dict)
    fillheadprint = defaultdict(dict)

    ioridict = defaultdict( lambda:defaultdict(dict))

    realneardict = gather_RSS_proximity(sources, "Unique_RSS_Site_files")

    #First we annotate each individual break in a per-patient manner, producing a folder full of individual DNA Damage files.
    # Each INDIVIDUAL SV will be annotated with the likely mutational source 
    for patientfile in glob.glob(f'Individual_Patient_SVs_ByRSSstatus/*_allSVs.csv'):
        patientidentifier = patientfile.split('/')
        patientidentifier = patientidentifier[1].split('_')
        filedatasource = patientidentifier[0]
        filedorr = patientidentifier[1]
        patientidentifier = patientidentifier[2]
        
        #this is necessary to separate diagnosis and relapse samples from one patient,
        #as well as samples for one patient from multiple data sources
        bigidentifier = patientidentifier + '_' + filedorr + '_' + filedatasource

        if bigidentifier not in individualpatientdict:
            individualpatientdict[bigidentifier] = []


        with open(patientfile) as infile, open(f'FullSVs_DNADamageExpanded/SplitByRSSProximity/{bigidentifier}_DNADamageExpanded.csv', 'w+') as individualoutfile:
            patientreader = csv.reader(infile)
            OUTbypatient = csv.writer(individualoutfile)

            OUTbypatientheader = ['Chromosome 1', 'Breakpoint 1', 'Chromosome 2', 'Breakpoint 2', 'Proximity to RSS' ,' ', 'Inter- or Intra-Chromosomal', 'Structural Variant Type', ' ', 'Chromosome of putative insertion', 'Site of Putative Insertion', 'Insertion site near an RSS?', 'Inserted Chromosome', 'Inserted Sequence Start' , 'Inserted Sequence End', 'Insertion Length', 'Relevant loci near insertion (LEFT)', 'Relevant loci near insertion (RIGHT)', ' ', 'Deletion Length']
            OUTbypatient.writerow(OUTbypatientheader)

            i = 0
            svdict = {}
            insertsites = []

            for line in patientreader:

                #Ensures consistent formatting of chromosome numbers
                chrom1 = str(line[0])
                if "chr" not in chrom1:
                    chrom1 = "chr" + chrom1
                pos1 = int(line[1])
                chrom2 = str(line[2])
                if "chr" not in chrom2:
                    chrom2 = "chr" + chrom2
                pos2 = int(line[3])
                proximity = line[4]

                #defines whether each SV is inter- or intra-chromosomal
                intra_or_inter = 'Interchromosomal'
                if chrom1 == chrom2:
                    intra_or_inter = 'Intrachromosomal'
                
                type = ''
                
                if i > 0:
                    #this is inefficient due to repeatedly iterating over the full list of SVs per pateint, 
                    # but for the scale of our datasets (approx 600 individual cases) still executes in under 30s
                    #This code is necessary as putative fragments of an insertion may be separated throughout the file
                    for entry in range(i):
                        type = ''
                        origchrom1 = str(svdict[entry][0])
                        origpos1 = int(svdict[entry][1])
                        origchrom2 = str(svdict[entry][2])
                        origpos2 = int(svdict[entry][3])

                        matched = False

                        #This next block of code defines insertions when two seperate SVs have breakpoints within
                        #50bp of one another, with a chunk of foreign DNA between.

                        #When inserted DNA appears to be from multiple chromosomes, sequences are instead defined 
                        #as Complex Insertions - I do not believe the exact mechanisms of these can be unpicked
                        #further based purely on chromosomal coordinates included in TARGET datasets.
                        if (pos1 - 50) <= origpos1 <= (pos1 + 50) and chrom1 == origchrom1:
                            insertsite = f'{pos1} - {origpos1}'
                            intochr = chrom1

                            if chrom2 != origchrom2:
                                type = 'Complex insertion'
                                insertsites.append(insertsite)
                                pass

                            else:
                                insertionfrom = chrom2
                                insertionstart = origpos2
                                insertionend = pos2
                                matched = True

                        elif (pos1 - 50) <= origpos2 <= (pos1 + 50) and chrom1 == origchrom2:
                            insertsite = f'{pos1} - {origpos2}'
                            intochr = chrom1

                            if chrom2 != origchrom1:
                                type = 'Complex insertion'
                                insertsites.append(insertsite)
                                pass

                            else:
                                insertionfrom = chrom1
                                insertionstart = origpos1
                                insertionend = pos2
                                matched = True

                        elif (pos2 - 50) <= origpos1 <= (pos2 + 50) and chrom2 == origchrom1:
                            insertsite = f'{pos2} - {origpos1}'
                            intochr = chrom2

                            if chrom1 != origchrom2:
                                type = 'Complex insertion'
                                insertsites.append(insertsite)
                                pass

                            else:
                                insertionfrom = chrom2
                                insertionstart = origpos2
                                insertionend = pos1
                                matched = True

                        elif (pos2 - 50) <= origpos2 <= (pos2 + 50) and chrom2 == origchrom2:
                            insertsite = f'{pos2} - {origpos2}'
                            intochr = chrom2

                            if chrom1 != origchrom1:
                                type = 'Complex insertion'
                                insertsites.append(insertsite)
                                pass

                            else:
                                insertionfrom = chrom1
                                insertionstart = origpos1
                                insertionend = pos1
                                matched = True

                        #If two SVs appear to constitute an insertion, we perform some further analysis
                        if matched == True:
                            insertionlength = abs(insertionstart - insertionend)

                            #if the insertion length appears to be under 100 bp, I discount the insertion - this is likely
                            #instead due to misalignment of a deletion or translocation

                            #Note that this is also likely to filter out inversions, but it would be unfeasible to
                            #attempt to define these seperately
                            if insertionlength <= 100:
                                break
                            else:
                                type = 'Insertion'

                            #runs locifinder function on the start and end position of the insertion
                            startloci = locifind(insertionstart, insertionfrom)
                            endloci = locifind(insertionend, insertionfrom)
                            bothloci = [startloci,endloci]

                            #I have elected to be liberal here in defining whether a sequence derives from Ig/Tcr
                            #If EITHER insertion end maps to this region, the insertion is defined as Ig/Tcr
                            if startloci != endloci:
                                fullloci = min(bothloci, key=len)
                            else:
                                fullloci = startloci    

                            #counts total SV loci in each patient
                            if fullloci not in locicountdict[filedatasource]:
                                locicountdict[filedatasource][fullloci] = 1
                            else:
                                locicountdict[filedatasource][fullloci] += 1  

                            #counts intra- and inter-chromosomal events in each patient
                            for iori in chromosomalchar:
                                if iori == intra_or_inter:
                                    if fullloci not in ioridict[filedatasource][iori]:
                                        ioridict[filedatasource][iori][fullloci] = 1
                                    else:
                                        ioridict[filedatasource][iori][fullloci] += 1

                            #Outputs all relevant inserts from a loci of interest to its own file (no header).
                            #This is needed for downstream loci_scanner program
                            if fullloci != 'No match to Ig or Tcr regions':
                                with open(f'FullSVs_DNADamageExpanded/Inserts_By_Loci/{fullloci}_{filedatasource}.csv', 'a+') as insertsout:
                                    OUTspecificloci = csv.writer(insertsout)
                                    specificlocirow = [bigidentifier, insertionfrom, insertionstart,insertionend, startloci, endloci]

                                    OUTspecificloci.writerow(specificlocirow)

                            insertsitessplit = insertsite.split(' - ')

                            #Determines whether the insertion site is near an RSS
                            try:
                                if realneardict[filedatasource][int(insertsitessplit[0])] == intochr:
                                    insertionsitenearrss = 'Yes'
                            except:
                                try:
                                    if realneardict[filedatasource][int(insertsitessplit[1])] == intochr:
                                        insertionsitenearrss = 'Yes'
                                except:
                                    insertionsitenearrss = 'No'

                            #Generates a BIG csv with intentional gaps - these are filled in by the reintegration_check program
                            with open(f'FullSVs_DNADamageExpanded/Inserts_By_Loci/INSERTION_INFORMATION_SUMMARY_{fullloci}_{filedatasource}.csv', 'a+') as fillableout:
                                OUTfillable = csv.writer(fillableout)       
                                fillablehead = ['Patient ID', 'Insertion site Chrom', 'Insertion Site', 'Insertion near an RSS?', 'Insertion near a gene segment?', 'Which gene segment?', 'Intragenic?', ' ', 'Insertion from Chrom', 'Insertion from loci', 'Insertion start site', 'Insertion start near RSS?','Insertion RSS in right orientation?', 'Insertion Start near gene segment?', 'Which gene segment?', 'Insertion Start in right orientation to gene seg?', ' ', 'Insertion end site', 'Insertion end near an RSS?', 'Insertion end RSS in right orientation?', 'Insertion end near gene segment?', ' Which Gene Segment?', ' Insertion end in right orientation to gene seg?', '  ', 'Both insert ends near RSSs?', 'Both insert ends near gene segments?', 'Putative Reintegration?']

                                #ensures the header line is only printed once
                                if fullloci not in fillheadprint[filedatasource]:
                                    OUTfillable.writerow(fillablehead)          
                                    fillheadprint[filedatasource][fullloci] = 1
                                
                                #checks whether the start of the inserted sequence is near an RSS
                                try:
                                    if realneardict[filedatasource][int(insertionstart)] == insertionfrom:
                                        instartnear = 'Yes'
                                    else: 
                                        instartnear = 'No'
                                except:
                                        instartnear = 'No'
                                
                                #checks whether the end of the inserted sequence is near an RSS
                                try:
                                    if realneardict[filedatasource][int(insertionend)] == insertionfrom:
                                        inendnear = 'Yes'
                                    else: 
                                        inendnear = 'No'
                                except:
                                        inendnear = 'No'

                                #checks whether BOTH insert ends are near RSSs 
                                if instartnear == 'Yes' and inendnear == 'Yes':
                                    bothendsnear = 'Yes'
                                else:
                                    bothendsnear = 'No'

                                filloutrow = [bigidentifier, intochr, insertsite, insertionsitenearrss, '  ','  ', '   ', '   ', insertionfrom, fullloci, insertionstart, instartnear, '   ','   ', '   ', '   ', ' ', insertionend, inendnear, '   ','   ', '   ', '  ', ' ', bothendsnear, '   ' , '   ']
                                OUTfillable.writerow(filloutrow)

                            extendit = [type,' ',intochr, insertsite, insertionsitenearrss, insertionfrom, insertionstart, insertionend, insertionlength, startloci, endloci]

                            #This code accounts for cases in which an SV corresponds to multiple insertions
                            #Allowing both to be annotated seperately. THIS IS DIFFERENT FROM A COMPLEX INSERTION
                            #Though may suggest presence of multiple clonal origins within samples.
                            if len(svdict[entry]) < 9:
                                svdict[entry].extend(extendit)
                            else:

                                svdict[entry][9] = str(svdict[entry][9]) + f'/{intochr}'
                                svdict[entry][10] = str(svdict[entry][10]) + f'/{insertsite}'
                                svdict[entry][11] = str(svdict[entry][11]) + f'/{insertionsitenearrss}'
                                svdict[entry][12] = str(svdict[entry][12]) + f'/{insertionfrom}'
                                svdict[entry][13] = str(svdict[entry][13]) + f'/{insertionstart}'
                                svdict[entry][14] = str(svdict[entry][14]) + f'/{insertionend}'
                                svdict[entry][15] = str(svdict[entry][15]) + f'/{insertionlength}'
                                svdict[entry][16] = str(svdict[entry][16]) + f'/{startloci}'
                                svdict[entry][17] = str(svdict[entry][17]) + f'/{endloci}'
                            break
                        
                        #Naturally, less info can be gathered for a Complex Insertion - No way to define start,
                        #end or length.
                        if type == 'Complex insertion':
                            extendit = [type,' ',intochr, insertsite, insertionsitenearrss, '-', '-', '-', '-', '-', '-']

                            if len(svdict[entry]) < 9:
                                svdict[entry].extend(extendit)
                            else:
                                svdict[entry][9] = str(svdict[entry][9]) + f'/{intochr}'
                                svdict[entry][10] = str(svdict[entry][10]) + f'/{insertsite}'
                                svdict[entry][11] = str(svdict[entry][11]) + f'/{insertionsitenearrss}'
                            break
                
                #Each patient SV added to a dictionary which is then iterated over further
                appendingline = [chrom1, pos1, chrom2, pos2, proximity, ' ', intra_or_inter]       
                if type == 'Insertion' or type == 'Complex insertion':
                    appendingline.extend(extendit)

                svdict[i] = appendingline
                i += 1
            
            #This 'n' variable is used for repeated iteration over the svdict dictionary
            #This ensures that anything contributing to a Complex Insertion can be appropriately annotated
            #Even if the "Complex" part occurs lower in the file
            n = 0

            while n <= 1:
                for entry in svdict:
                    try:
                        check = False
                        insertsitefromdict =  svdict[entry][10]
                        
                        if '/' in insertsitefromdict:
                            fromdictsplit = insertsitefromdict.split('/')

                            for x in fromdictsplit:
                                if x in insertsites:
                                    svdict[entry][7] = 'Complex insertion'
                                    svdict[entry][12], svdict[entry][13], svdict[entry][14], svdict[entry][15],svdict[entry][16], svdict[entry][17] = '-', '-', '-', '-', '-', '-'
                                    check = True
                            
                            if check == True: 
                                for y in fromdictsplit:
                                    if y not in insertsites:
                                        insertsites.append(y)
                        
                        else:
                            if insertsitefromdict in insertsites:
                                svdict[entry][7] = 'Complex insertion'
                                svdict[entry][12], svdict[entry][13], svdict[entry][14], svdict[entry][15], svdict[entry][16], svdict[entry][17] = '-', '-', '-', '-', '-', '-'
                    
                    except:
                        pass

                n += 1

            #For each individual SV, we ensure appropriate annotation with the type of mutation
            #and fill the remaining columns. These are output in individual files for each level of RSS proximity
            #along with a 'Compiled file' containing ALL SVs for each patient
            for entry in range(i):
                #Entries that were previously annotated as 'Insertion' or 'Complex insertion' are blanked for 
                #Deletion characteristics
                if len(svdict[entry]) > 9:
                    extendnodel = [' ', '-']
                    svdict[entry].extend(extendnodel)
                    OUTbypatient.writerow(svdict[entry])

                    tocompilefile = []  
                    tocompilefile.extend(svdict[entry])
                    individualpatientdict[bigidentifier].append(tocompilefile)

                else:
                    #Remaining entries which are intrachromosomal are defined as deletions.
                    #Deletion lengths (absolute distances between breakpoints) are recorded
                    if svdict[entry][6] == 'Intrachromosomal':
                        dellength = abs(svdict[entry][1] - svdict[entry][3])
                        extendnodel = [' ', dellength]
                        svdict[entry].append('Deletion')

                    #Remaining SVs which comprise multiple chromosomes are annotated as Translocations
                    #and are blanked for all other characteristics - cannot derive a length etc
                    else:
                        svdict[entry].append('Translocation')
                        extendnodel = [' ', '-']              

                    #writes to individual distance file
                    extendnoinsert = [' ','-', '-','-', '-', '-', '-','-', '-', '-']
                    svdict[entry].extend(extendnoinsert)
                    svdict[entry].extend(extendnodel)
                    OUTbypatient.writerow(svdict[entry])

                    #appends to a dict for writing overall compiled patient file
                    tocompilefile = []
                    tocompilefile.extend(svdict[entry])
                    individualpatientdict[bigidentifier].append(tocompilefile)

                    

    #Writes compiled files for each individual patient
    for patient in individualpatientdict:
        with open(f'FullSVs_DNADamageExpanded/SplitByRSSProximity/ByPatientCompiled/{patient}.csv', 'w+') as patientcompiledfile:
            OUTcompiledpatient = csv.writer(patientcompiledfile)

            OUTbycompiledheader = ['Proximity to RSSs', 'Chromosome 1', 'Breakpoint 1', 'Chromosome 2', 'Breakpoint 2' ,' ', 'Inter- or Intra-Chromosomal', 'Structural Variant Type', ' ', 'Chromosome of putative insertion', 'Site of Putative Insertion', 'Insertion site near an RSS?', 'Inserted Chromosome', 'Inserted Sequence Start' , 'Inserted Sequence End', 'Insertion Length', 'Relevant loci near insertion', ' ', 'Deletion Length']
            OUTcompiledpatient.writerow(OUTbycompiledheader)

            for individualsv in individualpatientdict[patient]:
                OUTcompiledpatient.writerow(individualsv)

    #Writes a metadata file for each loci from each data source with info on the percentage of SV events which are
    # intra- vs inter-chromosomal 
    for datasource in sources:
        for iori in ioridict[datasource]:
            for locus in ioridict[datasource][iori]:
                with open(f'FullSVs_DNADamageExpanded/Inserts_By_Loci/{locus}_{datasource}_insertmetadata.txt', 'a+') as insertsout:
                    insertsout.write(f'{locicountdict[datasource][locus]} insertions from {locus} \n')
                    insertsout.write(f'{ioridict[datasource][iori][locus]} insertions are {iori} \n')
                    percentiori = percentage(ioridict[datasource][iori][locus], locicountdict[datasource][locus])
                    insertsout.write(f' {percentiori} % are {iori} \n \n')
                        
    #Next, we compile information from across patients from each data source.
    for datasource in sources:
        dictbypatient = {}

        #This will write a number of files, facilitating varied analysis. Each is described below
        with open(f'FullSVs_DNADamageExpanded/SplitByRSSProximity/Summaries/{datasource}_DNAdamage_summary.csv', 'w+') as summaryfile, open(f'FullSVs_DNADamageExpanded/SplitByRSSProximity/Summaries/{datasource}_Insertions.csv', 'w+') as insertfile, open(f'FullSVs_DNADamageExpanded//SplitByRSSProximity/Summaries/{datasource}_Deletions.csv', 'w+') as deletionfile, open(f'FullSVs_DNADamageExpanded/SplitByRSSProximity/Summaries/{datasource}_Translocations.csv', 'w+') as translocationfile, open(f'FullSVs_DNADamageExpanded/SplitByRSSProximity/Summaries/{datasource}_Complex_Insertions.csv', 'w+') as complexfile, open(f'FullSVs_DNADamageExpanded/SplitByRSSProximity/Summaries/Cut&RunAssessment_{datasource}.csv', 'w+') as cutandrunfile:
            #OUTsummary contains numbers for each category of mutation per patient PER LEVEL OF RSS PROXIMITY
            OUTsummary = csv.writer(summaryfile)
            summaryheader = ['Patient ID', 'Diagnosis or Relapse', ' ', 'Total SVs', 'Interchromosomal Events', 'Intrachromosomal Events',' ', 'Number of SVs explained by Insertions', '% SVs are insertions', 'Inserts near Ig/Tcr loci', '% Inserts near Ig/Tcr loci', 'Average Insertion Length', 'Median Insertion Length', 'Number of SVs explained by Putative Complex Insertions', '% SVs are Complex insertions', ' ', 'Number of Deletions', '% SVs are deletions', 'Average Deletion Length', 'Median Deletion Length', ' ', 'Number of Translocations', '% SVs are translocations']
            OUTsummary.writerow(summaryheader)
            
            #OUTinserts contains all information pertaining to each individual insertion - what this means is
            #BOTH SVs which suggest presence of an insertion (SVs 1 and 2 for simplicity)
            #From this, we can derive where each insertion takes place, along with noting when Ig/Tcr loci are inserted
            OUTinserts = csv.writer(insertfile)
            insertsheader = ['Patient ID', 'Diagnosis or Relapse',' ', 'SV1 Left Chr', 'SV1 Left Pos', 'SV1 Right Chr', 'SV1 Right Pos', 'SV1 Proximity to RSSs', ' ', 'SV2 Left Chr', 'SV2 Left Pos', 'SV2 Right Chr', 'SV2 Right Pos', 'SV2 Proximity to RSSs', ' ','Insertion into Chromosome', 'Putative Insertion Site', 'Insertion site near an RSS?', 'Putative Inserted Sequence Chromosome', 'Putative Inserted Sequence Start', 'Putative Inserted Sequence End', 'Putative Insertion Length', 'Relevant Loci in Insert (LEFT)','Relevant Loci in Insert (RIGHT)', ' ', 'CutandRun insertion', 'Insertion near loci of interest' ]
            OUTinserts.writerow(insertsheader)  

            #OUTdeletions contains all information pertaining to each individual deletion
            OUTdeletions = csv.writer(deletionfile)
            deletionsheader = ['Patient ID', 'Diagnosis or Relapse', 'Chromosome', 'Left Pos', 'Right Pos', 'Proximity to RSS', 'Putative Deletion Length']
            OUTdeletions.writerow(deletionsheader)

            #OUTtranslocations contains all information pertaining to each individual insertion 
            OUTtranslocations = csv.writer(translocationfile)
            translocationsheader = ['Patient ID', 'Diagnosis or Relapse', 'Left Chromosome', 'Left Pos', 'Right Chromosome', 'Right Pos', 'Proximity to RSS']
            OUTtranslocations.writerow(translocationsheader)

            #OUTcomplex contains all information pertaining to each individual COMPLEX insertion - as for standard
            #insertions this provides two SVs per putative insertion event
            OUTcomplex = csv.writer(complexfile)
            complexheader = ['Patient ID', 'Diagnosis or Relapse',' ', 'SV1 Left Chr', 'SV1 Left Pos', 'SV1 Right Chr', 'SV1 Right Pos', 'SV1 Proximity to RSS', ' ', 'SV2 Left Chr', 'SV2 Left Pos', 'SV2 Right Chr', 'SV2 Right Pos', 'SV2 Proximity to RSS', ' ', 'Insertion into Chromosome', 'Insertion near an RSS?', 'Putative Insertion Site', ]
            OUTcomplex.writerow(complexheader)

            #OUTcutandrun is a HUGE csv containing all relevant stats from OUTsummary, separated by the proximity
            #of a given break to cRSS sequences
            OUTcutandrun = csv.writer(cutandrunfile)
            cutandrunheader = ['Patient ID', 'Diagnosis or Relapse', 'Total SVs', 'Total SVs which are definitively not cut and run', 'Total SVs with breaks at one RSS', 'Total SVs with breaks at two RSSs', ' ', 'Insertions not at RSSs', 'Percent total Svs are insertions not at RSSs', 'Percentage SVs not at RSSs which are insertions', 'Insertions not at RSSs with Ig/Tcr inserts', 'Percentage insertions not at RSSs have Ig/Tcr inserts', ' ', 'Inserts at RSSs',  'Percent total Svs are insertions at an RSS', 'Percent total SVs at one RSS which are insertions', 'Insertions at RSSs with Ig/Tcr inserts', 'Percentage at RSSs have Ig/Tcr inserts', ' ', 
                                'Deletions not near an RSS', 'Percentage of total SVs are deletions not near RSS', 'Percentage of SVs not near RSSs which are deletions', 'Deletions near one RSS', 'Percentage of total SVs which are deletions near one RSS', 'Percentage of SVs near one RSS which are deletions', 'Deletions near two RSSs', 'Percentage of total SVs are deletions near two RSSs', 'Percentage of SVs near two RSSs which are deletions', ' ',
                                'Translocations not near an RSS', 'Percentage of total SVs are translocations not near RSS', 'Percentage of SVs not near an RSS are translocations', 'Translocations near one RSS', 'Percentage of total SVs which are translocations near one RSS', 'Percentage of SVs near one RSS which are translocations', 'Translocations near two RSSs', 'Percentage of total SVs are translocations near two RSSs', 'Percentage of SVs near two RSSs which are translocations', ' ',
                                'Complex Insertions not near an RSS', 'Percentage of total SVs are Complex Insertions not near RSS', 'Percentage of SVs not near an RSS are Complex Insertions', 'Complex Insertions near one RSS', 'Percentage of total SVs which are Complex Insertions near one RSS', 'Percentage of SVs near one RSS which are Complex Insertions']
            OUTcutandrun.writerow(cutandrunheader)

            #We iterate over previously generated files to populate these files
            for damagefile in glob.glob(f'FullSVs_DNADamageExpanded/SplitByRSSProximity/*{datasource}_DNADamageExpanded.csv'):
                
                #Insertions are quantified as floats to account for the fact that each insertion event requires two SVs
                totalinsertions = 0.0
                totaldeletions = 0
                totaltranslocations = 0
                totalsvs = 0
                totalcomplex = 0
                insertsnearrelevantseqs = 0.0

                intrachromosomal = 0
                interchromosomal = 0

                #These variables hold all the information for various stats at various distances from RSSs
                totalsvs_0, totalsvs_1, totalsvs_2 = 0, 0, 0
                insertions_0, insertions_1 = 0, 0
                insertionsnear_0, insertionsnear_1 = 0, 0
                deletionsnear_0, deletionsnear_1, deletionsnear_2 = 0, 0, 0
                translocationsnear_0, translocationsnear_1, translocationsnear_2 = 0, 0, 0
                complex_0, complex_1 = 0, 0

                deletionlengths = []
                insertlengths = []

                patientidentifierlong = damagefile.split('/')
                patientidentifierlong = patientidentifierlong[2].split('_')
                patientidentifier = patientidentifierlong[0]
                patientdorr = patientidentifierlong[1]

                bigidentifier = patientidentifier + '_' + patientdorr

                insertinfodict = {}
                fragsdict = {}
                posdict = {}
                chrombypos = {}

                if bigidentifier not in dictbypatient:
                    dictbypatient[bigidentifier] = {}
            
                with open(damagefile) as patientdamage:

                    damagereader = csv.reader(patientdamage)
                    next(damagereader)

                    #iterates over each SV in the patient, incrementing various counters to generate final interesting numbers
                    for line in damagereader:
                        totalsvs += 1

                        iori = line[6]
                        if iori == 'Interchromosomal':
                            interchromosomal += 1
                        else:
                            intrachromosomal += 1
                        
                        sv1leftchr = line[0]
                        sv1leftpos = line[1]
                        sv1rightchr = line[2]
                        sv1rightpos = line[3]
                        sv1proximity = line[4]

                        damagetype = line[7]

                        if damagetype == 'Insertion':
                            linefrag = [sv1leftchr, sv1leftpos, sv1rightchr, sv1rightpos, sv1proximity, ' ']
                            totalinsertions += 1

                            #This is necessary to parse SVs involved in multiple insertion events.
                            if '/' in line[9]:
                                insertedinto = line[9].split('/')
                                insertedpos = line[10].split('/')
                                insertnearrss = line[11].split('/')
                                insertedchrs = line[12].split('/')
                                insertedstart = line[13].split('/')
                                insertedend = line[14].split('/')
                                insertedlengths = line[15].split('/')
                                insertedloci1 = line[16].split('/')
                                insertedloci2 = line[17].split('/')
                            
                                for length in insertedlengths:
                                    if length not in insertlengths:
                                        insertlengths.append(int(length))

                                for x in range(len(insertedchrs)):
                                    tempinsertinto = insertedinto[x]
                                    tempinsertpos = insertedpos[x]
                                    tempinsertnearrss = insertnearrss[x]
                                    tempinsertchr = insertedchrs[x]
                                    tempinsertstart = insertedstart[x]
                                    tempinsertend = insertedend[x]
                                    tempinsertlen = int(insertedlengths[x])
                                    tempinsertlocus1 = insertedloci1[x]
                                    tempinsertlocus2 = insertedloci2[x]

                                    if tempinsertlocus1 != 'No match to Ig or Tcr regions' or tempinsertlocus2 != 'No match to Ig or Tcr regions' :
                                        insertsnearrelevantseqs += 0.5
                                        nearloi = 1
                                    else:
                                        nearloi = 0
                                    
                                    #if an insertion site is near an RSS, we classify this as "putative cut-and-run",
                                    #although naturally other sources of damage may be involved.
                                    if tempinsertnearrss == 'Yes':
                                        putativecutandrun = 1
                                        totalsvs_1 += 1
                                        insertions_1  += 1
                                        if nearloi == 1:
                                            insertionsnear_1 += 1
                                    else:
                                        putativecutandrun = 0
                                        totalsvs_0 += 1
                                        insertions_0 += 1
                                        if nearloi == 1:
                                            insertionsnear_0 += 1
                                    
                                    #the length of the insertion is used as a key for the dictionary of line fragments
                                    #This has yet to cause a problem, as we have checked these are unique through out datasets
                                    #Though it is feasible that this COULD cause issues in larger datasets.
                                    if tempinsertlen not in fragsdict:
                                        fragsdict[tempinsertlen] = []
                                    if linefrag not in fragsdict[tempinsertlen]:
                                        fragsdict[tempinsertlen].append(linefrag)

                                    insertinfo = [tempinsertinto, tempinsertpos, tempinsertnearrss, tempinsertchr, tempinsertstart, tempinsertend, tempinsertlen, tempinsertlocus1, tempinsertlocus2, ' ', putativecutandrun, nearloi]
                                    insertinfodict[tempinsertlen] = insertinfo
                            
                            else:
                                #same process applied for SV breakpoints only involved in one event
                                if line[15] not in insertlengths:
                                    insertlengths.append(int(line[15]))
                                
                                insertioninto = line[9]
                                insertedpos = line[10]
                                insertionnearrss = line[11]
                                insertedchr = line[12]
                                insertedstart = line[13]
                                insertedend = line[14]
                                insertedlength = int(line[15])
                                insertedlocus1 = line[16]
                                insertedlocus2 = line[17]

                                if insertedlength not in fragsdict:
                                        fragsdict[insertedlength] = []
                                
                                if insertedlocus1 != 'No match to Ig or Tcr regions' or insertedlocus2 != 'No match to Ig or Tcr regions':
                                    insertsnearrelevantseqs += 1
                                    nearloi = 1
                                else:
                                    nearloi = 0

                                if insertionnearrss == 'Yes':
                                    putativecutandrun = 1
                                    totalsvs_1 += 1
                                    insertions_1  += 1
                                    if nearloi == 1:
                                        insertionsnear_1 += 1
                                else:
                                    putativecutandrun = 0
                                    totalsvs_0 += 1
                                    insertions_0 += 1
                                    if nearloi == 1:
                                        insertionsnear_0 += 1

                                if linefrag not in fragsdict[insertedlength]:
                                    fragsdict[insertedlength].append(linefrag)

                                insertinfo = [insertioninto,insertedpos, insertionnearrss, insertedchr,insertedstart, insertedend, insertedlength, insertedlocus1, insertedlocus2, ' ', putativecutandrun, nearloi]
                                insertinfodict[insertedlength] = insertinfo

                        elif damagetype == 'Deletion':
                            totaldeletions += 1
                            deletionlengths.append(int(line[19]))

                            #Writes each event to the Deletion file
                            deloutline = [patientidentifier, patientdorr, sv1leftchr, sv1leftpos, sv1rightpos, sv1proximity, line[19]]
                            OUTdeletions.writerow(deloutline)

                            #Deletions separated by proximity of breakpoints to RSSs
                            if 'Neither' in sv1proximity:
                                totalsvs_0 += 1
                                deletionsnear_0 += 1
                            elif 'One' in sv1proximity:
                                totalsvs_1 += 1
                                deletionsnear_1 += 1
                            elif 'Both' in sv1proximity:
                                totalsvs_2 += 1
                                deletionsnear_2 += 1
                        
                        elif damagetype == 'Translocation':
                            totaltranslocations += 1

                            #Writes each event to the Translocation file
                            transoutline = [patientidentifier, patientdorr, sv1leftchr, sv1leftpos, sv1rightchr, sv1rightpos, sv1proximity]
                            OUTtranslocations.writerow(transoutline)
                            
                            #translocations separated by proximity of breakpoints to RSSs
                            if 'Neither' in sv1proximity:
                                totalsvs_0 += 1
                                translocationsnear_0 += 1
                            elif 'One' in sv1proximity:
                                totalsvs_1 += 1
                                translocationsnear_1 += 1
                            elif 'Both' in sv1proximity:
                                totalsvs_2 += 1
                                translocationsnear_2 += 1
                        
                        elif damagetype == 'Complex insertion':
                            totalcomplex += 1
                            linefrag = [sv1leftchr, int(sv1leftpos), sv1rightchr, int(sv1rightpos), sv1proximity, ' ']

                            #This is necessary to parse SVs involved in multiple COMPLEX insertion events.
                            if '/' in line[8]:
                                insertedinto = line[9].split('/')
                                insertedpos = line[10].split('/')
                                insertnearrss = line[11].split('/')

                                for x in range(len(insertedinto)):
                                    tempinsertinto = insertedinto[x]
                                    tempinsertpos = insertedpos[x]
                                    tempinsertnearrss = insertnearrss[x]
                                    
                                    if tempinsertpos not in posdict:
                                        posdict[tempinsertpos] = []

                                    #As for normal insertions, complex insertions can only be near or not near
                                    #an RSS at the insertion site
                                    if tempinsertnearrss == 'Yes':
                                        totalsvs_1 += 1
                                        complex_1  += 1
                                    else:
                                        totalsvs_0 += 1
                                        complex_0 += 1

                                    if tempinsertpos not in chrombypos:
                                        chrombypos[tempinsertpos] = [tempinsertinto,tempinsertnearrss]

                                    posdict[tempinsertpos].append(linefrag)
                            
                            else:
                                #single Complex Insertions handled in the same way
                                insertedinto = line[9]
                                insertedpos = line[10]
                                insertnearrss = line[11]

                                if insertedpos not in posdict:
                                    posdict[insertedpos] = []
                                
                                if insertnearrss == 'Yes':
                                    totalsvs_1 += 1
                                    complex_1  += 1
                                else:
                                    totalsvs_0 += 1
                                    complex_0 += 1

                                if insertedpos not in chrombypos:
                                    chrombypos[insertedpos] = [insertedinto, insertnearrss]

                                posdict[insertedpos].append(linefrag)
                    
                    #Assembles individual SV 'fragments' which compose an Insertion
                    for insertionlen in insertinfodict:
                        insertout = [patientidentifier, patientdorr, ' ']

                        for individualevent in fragsdict[insertionlen]:
                            insertout.extend(individualevent)
                        
                        insertout.extend(insertinfodict[insertionlen])
                        OUTinserts.writerow(insertout)
                    
                    #As before but for Complex Insertions. This takes a bit more work given the complex nature.
                    for insertionpos in posdict:
                        fraglist = []
                        complexout = [patientidentifier, patientdorr, ' ']

                        if len(posdict[insertionpos]) == 2:
                            for individualevent in posdict[insertionpos]:
                                complexout.extend(individualevent)

                            complexout.extend(chrombypos[insertionpos])
                            complexout.append(insertionpos)
                            OUTcomplex.writerow(complexout)

                        else:
                            #Here we need to iterate over all complex insertions AGAIN to make sure that every breakpoint is caught
                            for fragment in posdict[insertionpos]:
                                #Without this check that fragments aren't matching with themselves, lines were being assembled incorrectly
                                for prevfrag in fraglist:
                                    if fragment != prevfrag:
                                        check = False
                                        complexout = [patientidentifier, patientdorr, ' ']
                                        if fragment[0] == prevfrag[0] and (fragment[1] - 50) <= prevfrag[1] <= (fragment[1] + 50):
                                            check = True
                                        elif fragment[0] == prevfrag[2] and (fragment[1] - 50) <= prevfrag[3] <= (fragment[1] + 50):
                                            check = True
                                        elif fragment[2] == prevfrag[0] and (fragment[3] - 50) <= prevfrag[1] <= (fragment[3] + 50):
                                            check = True
                                        elif fragment[2] == prevfrag[2] and (fragment[3] - 50) <= prevfrag[3] <= (fragment[3] + 50):
                                            check = True
                                        
                                        if check == True:    
                                            complexout.extend(fragment) 
                                            complexout.extend(prevfrag)
                                            complexout.extend(chrombypos[insertionpos])
                                            complexout.append(insertionpos)
                                            OUTcomplex.writerow(complexout)

                                if fragment not in fraglist:        
                                    fraglist.append(fragment)
                    
                    #Some basic stats are calculated for insert lengths
                    try: 
                        meaninsertlen = statistics.mean(insertlengths)
                        medianinsertlen = statistics.median(insertlengths)
                    except:
                        meaninsertlen = 'N/A'
                        medianinsertlen = 'N/A'

                    #And for deletion lengths
                    try:
                        meandeletionlen = statistics.mean(deletionlengths)
                        mediandeletionlen = statistics.median(deletionlengths)
                    except:
                        meandeletionlen = 'N/A'
                        mediandeletionlen = 'N/A'

                    #We calculate the percentage of total SVs at a given level of RSS proximity which correspond to
                    #each type of mutation
                    if totalsvs != 0:
                        percentareinsertions = percentage(totalinsertions,totalsvs)
                        percentarecomplex = percentage(totalcomplex,totalsvs)
                        percentaredels = percentage(totaldeletions,totalsvs)
                        percentaretranslocations = percentage(totaltranslocations,totalsvs)
                    else:
                        percentareinsertions = 'N/A'
                        percentarecomplex = 'N/A'
                        percentaredels = 'N/A'
                        percentaretranslocations = 'N/A'
                    
                    #Exhaustive block to calculate percentages near NO RSSs
                    if totalsvs_0 != 0:
                        percenttotalinsert_0 = percentage(insertions_0,totalsvs)
                        percentinsert_0 = percentage(insertions_0,totalsvs_0)
                        percenttotaldeletions_0 = percentage(deletionsnear_0, totalsvs)
                        percentdeletionsnear_0 = percentage(deletionsnear_0, totalsvs_0)
                        percenttotaltranslocations_0 = percentage(translocationsnear_0,totalsvs)
                        percenttranslocationsnear_0 = percentage(translocationsnear_0,totalsvs_0)
                        percenttotalcomplex_0 = percentage(complex_0,totalsvs)
                        percentcomplex_0 = percentage(complex_0,totalsvs_0)
                        if insertions_0 != 0:
                            percentnearloci_0 = percentage(insertionsnear_0, insertions_0)
                        else:
                            percentnearloci_0 = 'N/A'
                    else:
                        percenttotalinsert_0 = 'N/A'
                        percentinsert_0 = 'N/A'
                        percenttotaldeletions_0 = 'N/A'
                        percentdeletionsnear_0 = 'N/A'
                        percenttotaltranslocations_0 = 'N/A'
                        percenttranslocationsnear_0 = 'N/A'
                        percenttotalcomplex_0 = 'N/A'
                        percentcomplex_0 = 'N/A'
                        percentnearloci_0 = 'N/A'

                    #Exhaustive block to calculate percentages near ONE RSS
                    if totalsvs_1 != 0:
                        percenttotalinsert_1 = percentage(insertions_1,totalsvs)
                        percentinsert_1 = percentage(insertions_1,totalsvs_1)
                        percenttotaldeletions_1 = percentage(deletionsnear_1, totalsvs)
                        percentdeletionsnear_1 = percentage(deletionsnear_1, totalsvs_1)
                        percenttotaltranslocations_1 = percentage(translocationsnear_1,totalsvs)
                        percenttranslocationsnear_1 = percentage(translocationsnear_1,totalsvs_1)
                        percenttotalcomplex_1 = percentage(complex_1,totalsvs)
                        percentcomplex_1 = percentage(complex_1,totalsvs_1)
                        if insertions_1 != 0:
                            percentnearloci_1 = percentage(insertionsnear_1, insertions_1)
                        else:
                            percentnearloci_1 = 'N/A'
                    else:
                        percenttotalinsert_1 = 'N/A'
                        percentinsert_1 = 'N/A'
                        percenttotaldeletions_1= 'N/A'
                        percentdeletionsnear_1 = 'N/A'
                        percenttotaltranslocations_1 = 'N/A'
                        percenttranslocationsnear_1 = 'N/A'
                        percenttotalcomplex_1 = 'N/A'
                        percentcomplex_1 = 'N/A'
                        percentnearloci_1 = 'N/A'
                    
                    #Exhaustive block to calculate percentages near TWO RSSs
                    if totalsvs_2 != 0:
                        percenttotaldeletions_2 = percentage(deletionsnear_2, totalsvs)
                        percentdeletionsnear_2 = percentage(deletionsnear_2, totalsvs_2)
                        percenttotaltranslocations_2 = percentage(translocationsnear_2,totalsvs)
                        percenttranslocationsnear_2 = percentage(translocationsnear_2,totalsvs_2)
                    else:
                        percenttotaldeletions_2= 'N/A'
                        percentdeletionsnear_2 = 'N/A'
                        percenttotaltranslocations_2 = 'N/A'
                        percenttranslocationsnear_2 = 'N/A'

                    #We also calculate the percentage of insertions which are near our loci of interest
                    try:
                        percentnearloci = (insertsnearrelevantseqs/totalinsertions)*100
                    except:
                        percentnearloci = 'N/A'

                    wholesummaryoutrow = [patientidentifier, patientdorr,' ', totalsvs, interchromosomal, intrachromosomal, ' ', totalinsertions, percentareinsertions, insertsnearrelevantseqs, percentnearloci, meaninsertlen, medianinsertlen, totalcomplex, percentarecomplex, ' ', totaldeletions, percentaredels, meandeletionlen, mediandeletionlen, ' ', totaltranslocations, percentaretranslocations]
                    OUTsummary.writerow(wholesummaryoutrow)

                    #Constructs and prints the output to the Cut-and-Run file.
                    cutandrunoutrow = [patientidentifier, patientdorr, totalsvs, totalsvs_0, totalsvs_1, totalsvs_2, ' ', insertions_0, percenttotalinsert_0, percentinsert_0, insertionsnear_0, percentnearloci_0, ' ', 
                                        insertions_1, percenttotalinsert_1, percentinsert_1, insertionsnear_1, percentnearloci_1, ' ',
                                        deletionsnear_0, percenttotaldeletions_0, percentdeletionsnear_0, deletionsnear_1, percenttotaldeletions_1, percentdeletionsnear_1, deletionsnear_2, percenttotaldeletions_2, percentdeletionsnear_2,' ',
                                        translocationsnear_0, percenttotaltranslocations_0, percenttranslocationsnear_0, translocationsnear_1,percenttotaltranslocations_1, percenttranslocationsnear_1, translocationsnear_2, percenttotaltranslocations_2, percenttranslocationsnear_2, ' ',
                                        complex_0, percenttotalcomplex_0, percentcomplex_0, complex_1, percenttotalcomplex_1, percentcomplex_1 ]
                    OUTcutandrun.writerow(cutandrunoutrow)


main()