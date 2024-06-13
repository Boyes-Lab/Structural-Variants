import argparse
import csv
import glob
import os
import psutil
import subprocess as sub
import time
import urllib.request
import zipfile

from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.ui import Select
from selenium.webdriver.common.by import By

#global variable to count number of translocation events
total_count = 0
rand_count = 0

#non-essential function to prepare folders needed for running, for portability.
def folderprep():
    folders_needed = ['Bed_files', 'Random_beds', 'Fasta_files', 'Random_fastas', 'RSS_files', 'Random_rss', 'Unique_RSS_Site_files', 'Random_unique_sites', 'Mutations_near_RSSs', 'Lengths_of_Deletions','Randoms_near_RSSs', 'OUTPUTS']
    for folder_name in folders_needed:
        try:
            os.mkdir(folder_name)
        except:
            pass
    
    for file in glob.glob('Bed_files/*'):
        print(file)
        os.remove(file)

class Breakfile():

    def __init__(self, file_name):

        self.file_name = file_name
        self.name = os.path.splitext(os.path.basename(self.file_name))[0]
        
        self.neo_bed_filename = f'Bed_files/{self.name}.bed'
        self.fasta_file = f'Fasta_files/{self.name}.fasta'
        self.rss_file = f'RSS_files/{self.name}.rss.txt'
        self.unique_rss = f'Unique_RSS_Site_files/{self.name}.csv'
        self.mutations = f'Mutations_near_RSSs/{self.name}.csv'
        self.length = f'Lengths_of_Deletions/{self.name}.csv'

    #converts input files into workflow-compliant bed files and uses each breakpoint to generate nearby "random" sites for comparison
    def make_bed(self, window, random, orig_done):
        bed_file = open(self.neo_bed_filename, "a+")

        orig_bed_file = open(f'{self.name}.bed');
        filepartition = 1
        random_bed = open(f'Random_beds/{self.name}_{random}tobreak_RANDOM.bed', 'w+')
        linesinoutput = 0
        linesinrandoutput = 0
        outputs = [self.neo_bed_filename,]
        randoutputs = [f'Random_beds/{self.name}_{random}tobreak_RANDOM.bed',]

        for line in orig_bed_file:
            workingline = line.split('\t')
            chromosome = workingline[0]
            position = int(workingline[1])
            if orig_done == 0:
                bed_file.write('\t'.join([chromosome, str(position - window), str(position + window) +'\n']));
                linesinoutput += 1
                #Some of our input files were too large for RSS site to easily handle. This just breaks up input into managable 5000 sequence chunks
                if linesinoutput %5000 == 0:
                    filebreak = f'Bed_files/{self.name}_{filepartition}.bed'
                    outputs.append(filebreak)
                    bed_file = open(filebreak, 'w')
                global total_count
                total_count += 1 

            #writes a unique bed file containing sequences n bases downstream of target structural variant    
            random_bed.write('\t'.join([chromosome, str(position - window + random), str(position + window + random) +'\n']))
            linesinrandoutput += 1

            if linesinrandoutput % 5000 == 0:
                randfilebreak = f'Random_beds/{self.name}_{filepartition}_{random}tobreak_RANDOM.bed'
                randoutputs.append(randfilebreak)
                random_bed = open(randfilebreak, 'w')
                filepartition += 1

        orig_done += 1
        return outputs, randoutputs, orig_done

class RSSsite_interface():
    #Interfaces with the RSSsite to identify RSSs in .fasta files
    def __init__(self, fasta, output, rss_type):
        self.fasta = fasta
        self.output = output
        self.rss_type = {'12':'12 bases', '23': '23 bases', 'both': 'both'}.get(rss_type)
    
    def communicate(self):
        #Uploads the generated .fasta file to RSSsite and analyses the sequences for potential RSSs 

        options = Options()
        options.add_argument('-headless')
        upload = f'{os.getcwd()}/{self.fasta}'
        driver = webdriver.Firefox()
        driver.get('http://www.itb.cnr.it/rss/analyze.html')

        select_species_option = Select(driver.find_element("name","species"))
        select_species_option.select_by_visible_text('human')

        rss_choice = self.rss_type
        select_rss_option = Select(driver.find_element("name",'spacer'))
        select_rss_option.select_by_visible_text(str(rss_choice))

        upload_file = driver.find_element("name","upfile")
        upload_file.send_keys(upload)

        analyse_sequence_button = driver.find_element(By.XPATH,"/html/body/table[2]/tbody/tr/td[3]/table/tbody/tr[1]/td/table/tbody/tr[8]/td[2]/input")
        analyse_sequence_button.click()

        time.sleep(60)
        cpu = psutil.cpu_percent()
        while cpu > 20:
            cpu = psutil.cpu_percent()
        
        download_link = driver.find_element(By.LINK_TEXT, 'Click here to get the zipped txt tab separated version of this table')
        download_address = download_link.get_attribute('href')

        urllib.request.urlretrieve(download_address, self.output)
        driver.close()

        file_to_extract = zipfile.ZipFile(self.output)
        file_to_extract.extractall()
        files = glob.glob('file*.txt')
        rss_file = files.pop(0)
        os.rename(rss_file, self.output)

        for f in files:
            os.remove(f)
        
        return rss_file

def count_rss(rss_file, unique_rss):
    #Counts RSSs with a RIC score above the threshold
    
    def remove_fails(rss_file):
        with open(rss_file) as reader:
            rss = [x for x in reader if not 'FAIL' in x.split('\t')[-1]]
            rss = [x for x in rss if x[:3] == 'chr']  
        return rss

    rss_passes = remove_fails(rss_file)
    print(len(rss_passes))
    rss_unique = {x.split('\t')[0] for x in rss_passes}

    #Returns a 'Unique_RSS_Site' file listing the window used for detection of each putative cRSS
    header = ('Chromosome', 'Search Window Start', 'Search Window End')
    with open(unique_rss, 'w+', newline = '') as uniques_csv:
        writer = csv.writer(uniques_csv)
        writer.writerow(header)
        for unique in sorted(rss_unique):
            y = unique.replace('-', ':').split(':')
            writer.writerow(y)
        

    print('...unique RSS sites file made')
    return len(rss_unique)

def main():
    def get_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('-w','--window', default = 45, type = int, help = 'Defines the region around each breakpoint to be searched for RSSs, defailt of 45 bp')
        parser.add_argument('-RSS', '--RSS_type', default = 'both', help = 'Defines the type of RSS to search for at each breakpoint. Default is both 12 and 23 RSSs, but either may be specified alone')
        parser.add_argument('-m', '--mutations', default = 0, type = int, help = 'Option to find mutations near unique RSSs, enter 1 to run')
        parser.add_argument('-g', '--genome', help = 'Genome in fasta format')
        parser.add_argument('-rand', '--randomiser', default = (0,100,300,1000), type = list, help = 'number of basepairs downstream (+ve number) or upstream (-ve number) to take a so called random site')
        return parser.parse_args()

    args = get_args()
    case_no = 0

    with open ('OUTPUTS/Mutation_Count_Near_RSSs_By_Case.csv', mode = 'w+', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        firstrow = 0
        randlist = args.randomiser
        randfastalist = {}
        orig_done = 0

        #the inputs here are .bed files in cwd due to our processing pipeline for available data
        # other input types should be easily accessed with minimal modifications.
        for bed in sorted(glob.glob("*.bed")):
            global total_count
            total_count = 0
            global rand_count
            rand_count = 0

            case_no += 1
            print(f'Processing Case {case_no}')
            currentfile = Breakfile(bed)
            randbeddict = {}
            bedlist = []

            #for each designated distance, the equivalent pipeline will be run on sequences n basepairs downstream.
            for random in randlist:
                bedlistaddon, randbeddict[random], orig_done_update = currentfile.make_bed(args.window, random, orig_done)
                orig_done = orig_done_update
                for item in bedlistaddon:
                     if item not in bedlist:
                        bedlist.append(item)

            #sets up column names on output file
            header = ['Patient Identifier', 'Mutation near RSSs', 'Mutation not near RSSs', 'Total Mutations','% near RSSs', ' ']
            for random in randlist:
                additions = [f'Random sites {random} from breakpoints near RSSs', 'Random sites not near RSSs', 'Total Random sites', f'% Randoms {random} bp from breakpoints near RSSs', ' ']
                header.extend(additions)

            if firstrow == 0:
                csv_writer.writerow(header)
                firstrow += 1

            fastaout = currentfile.fasta_file
            partition = 1
            fastalist = []

            #makes one or more (based on length) fasta files of patient breakpoints
            for bedfile in bedlist:
                print(bedfile)
                cmd = f'bedtools getfasta -fi {args.genome} -bed {bedfile} -fo {fastaout} ' 
                sub.run(cmd.split())
                fastalist.append(fastaout)
                fastaout = f'Fasta_files/{currentfile.name}_{partition}.fasta'
                partition += 1

            print('...fasta files made')

            #makes one or more fasta files for EACH DISTANCE specified in randomisation files
            for random in randlist:
                randfastaout = f'Random_fastas/{currentfile.name}_{random}tobreak_RANDOM.fasta'
                partition = 1
                randfastalist[random] = []

                for randbed in randbeddict[random]:
                    cmd = f'bedtools getfasta -fi {args.genome} -bed {randbed} -fo {randfastaout} ' 
                    sub.run(cmd.split())
                    randfastalist[random].append(randfastaout)
                    randfastaout = f'Random_fastas/{currentfile.name}_{partition}_{random}tobreak_RANDOM.fasta'
                    partition += 1
                print(f'...{random} to random fasta file made')

            rssoutput = currentfile.rss_file
            uniqueoutput = currentfile.unique_rss
            filepartition = 1
            totalrsscount = 0

            #communicates with RSS site to find patient breakpoints near putative cRSSs, count these, and output sequences of each to a seperate file
            for fastafile in fastalist:
                
                interface = RSSsite_interface(fastafile, rssoutput, args.RSS_type)
                interface.communicate()
                rss_count = count_rss(rssoutput, uniqueoutput)
                totalrsscount += rss_count
                print(f'rsss counted for {fastafile} - {rss_count}')
                rssoutput = f'RSS_files/{currentfile.name}_{filepartition}.rss.txt'
                uniqueoutput = f'Unique_RSS_Site_files/{currentfile.name}_{filepartition}.csv'
                filepartition += 1

            total_rand_rsscount = {}

            #communicates with RSS site to find sites n basepairs from each patient breakpoint near putative cRSSs
            #these are counted, and sequences of each output to a seperate file
            for random in randlist:
                randrssoutput = f'Random_rss/{currentfile.name}_{random}tobreak_RANDOM.rss.txt'
                randuniqueoutput = f'Random_unique_sites/{currentfile.name}_{random}tobreak_RANDOM.csv'
                filepartition = 1


                for randfasta in randfastalist[random]:

                    interface = RSSsite_interface(randfasta, randrssoutput, args.RSS_type)
                    interface.communicate()
                    random_rss_count = count_rss(randrssoutput, randuniqueoutput)
                    try:
                        total_rand_rsscount[random] += random_rss_count
                    except:
                        total_rand_rsscount[random] = random_rss_count
                    print(f'rsss counted for {randfasta} - {random_rss_count}')
                    randrssoutput = f'Random_rss/{currentfile.name}_{filepartition}_{random}tobreak_RANDOM.rss.txt'
                    randuniqueoutput = f'Random_unique_sites/{currentfile.name}_{filepartition}_{random}tobreak_RANDOM.csv'
                    filepartition += 1


            #Outputs relevant figures to .csv file
            experiment = [f'{currentfile.name}', totalrsscount, (total_count - totalrsscount), total_count, ((totalrsscount/total_count)*100), ' ']
            for random in randlist:
                experimentaddons = [total_rand_rsscount[random], (total_count - total_rand_rsscount[random]), total_count,((total_rand_rsscount[random]/total_count)*100), ' ']
                experiment.extend(experimentaddons)
            csv_writer.writerow(experiment) 

start_time = time.time()    
folderprep()
main()
print('\n*RSS analysis complete*')

#this can be taken out if desired, just shows program runtime
print('--- %s seconds ---' % (time.time() - start_time))
        
