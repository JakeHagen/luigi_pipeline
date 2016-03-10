import luigi
import luigi.postgres
import psycopg2
import subprocess
import os                                         
import pandas as pd 
import datetime

class parameters(luigi.Config):
    ''' 
    Class to contain all parameters. Most of these will need to be set from luigi.cfg file in working directory, 
    which is good because it acts as a log for parameters of experiment
    '''
    
    fastqs = luigi.Parameter(default = None)
    exp_dir = luigi.Parameter(default = os.getcwd())
    genome_fasta = luigi.Parameter(default = None)
    genome_gtf = luigi.Parameter(default = None)
    star_genome_folder = luigi.Parameter(default = None)
    read_length = luigi.IntParameter(default = 100)
    stranded = luigi.BoolParameter(default = False)
    paried = luigi.IntParameter(default = 0)                                            
    cores = luigi.IntParameter(default = 6)    
    exp_name = luigi.Parameter(default = datetime.date.today().strftime("%B%d,%Y"))
    postgres_password = luigi.Parameter(default = None)
    star_genome_index = luigi.Parameter(default = None)    
    star_align = luigi.Parameter(default = None)
    star_command = luigi.Parameter(default = None)
    postgres_host = luigi.Parameter(default = 'localhost')
    postgres_database = luigi.Parameter(default = 'RNA')
    postgres_user = luigi.Parameter(default = 'hagenj02')
    

class fastqs(luigi.Task):
    '''Takes fastqs from parameters (specified in python.cfg) and
       returns dictionary with sample name and luigi.LocalTarget for the fastq
    '''

    def run(self):    
        fastq_dict = {}
        for line in parameters().fastqs.splitlines():
            sample, path = line.split(":")
            fastq_dict[sample] = path
        return fastq_dict
    
    def output(self):
        fastq_dict = self.run()
        return {sample:luigi.LocalTarget(path) for sample, path in fastq_dict.items()} 


class index_STAR_genome(luigi.Task):
    
    def run(self):
        if not os.path.exists(parameters().star_genome_folder):
            os.makedirs(parameters().star_genome_folder) 
        os.chdir(parameters().star_genome_folder)
        star_command = [
                        'STAR',
                        '--runThreadN %d' % parameters().cores,
                        '--runMode genomeGenerate',
                        '--genomeDir %s' % parameters().star_genome_folder,
                        '--genomeFastaFiles %s' % parameters().genome_fasta,
                        '--sjdbGTFfile %s' % parameters().genome_gtf,
                        '--sjdbOverhang %d' % (parameters().read_length - 1),
                        ]
        star = subprocess.Popen(star_command)
        star.wait()
        if star.returncode != 0:
            subprocess.call(['rm -rf %s' % parameters().star_genome_folder])
        os.rename('Log.out', '%s/Log.out' % parameters().star_genome_folder)
    
    def output(self):
        return luigi.LocalTarget('%s/Log.out' % parameters().star_genome_folder)


class exp_dir(luigi.Task):
    
    def requires(self):
        return index_STAR_genome()
    
    def run(self):
        os.makedirs(parameters().exp_dir)
    
    def output(self):
        return luigi.LocalTarget(parameters().exp_dir)


class star_align(luigi.Task):
    
    sample = luigi.Parameter()
    file_location = luigi.Parameter()

    def requires(self):
        return index_STAR_genome(), exp_dir(), fastqs()

    def run(self):
        if not os.path.exists('%s/%s/star' % (parameters().exp_dir, self.sample)):
            os.makedirs('%s/%s/star' % (parameters().exp_dir, self.sample))
        star_command = [
                        'STAR',
                        '--genomeDir %s' % parameters().star_genome_folder,
                        '--readFilesIn %s' % self.file_location.path,
                        '--runThreadN %d' % parameters().cores,
                        '--outFileNamePrefix %s/%s/star/%s.' %
                            (parameters().exp_dir, self.sample, self.sample)
                        ]
        #Append the extra parameters we wont from the config file
        for line in parameters().star_align.splitlines():
            star_command.append(line)
        star = subprocess.Popen(star_command)
        star.wait()

        #Below removes folder if star command was unsuccessful 
        #(files remain when star fails)
        #renaming bam insures that if batch job gets cut while star is running
        #the task will not be marked as complete 
        
        if star.returncode != 0:
            subprocess.call(['rm', '-rf', '%s/%s/star' 
                                % (parameters().exp_dir, self.sample)])
        else:
            os.rename('%s/%s/star/%s.Aligned.sortedByCoord.out.bam' 
                        % (parameters().exp_dir, self.sample, self.sample), 
                        '%s/%s/star/%s.bam' 
                            % (parameters().exp_dir, self.sample, self.sample))
    def output(self):
        return luigi.LocalTarget('%s/%s/star/%s.bam' 
                                    % (parameters().exp_dir, self.sample, self.sample))


class all_star_align(luigi.Task):
    '''Requires splits out samples and runs star align for eacg sample,
    output is dictionary with sample name and luigi target of bam file
    '''
    def requires(self):
        fastq_dict = fastqs().output()
        return {
                s:star_align(sample = s, file_location = p) 
                    for s,p in fastq_dict.items()
                }
    def output(self):       
        bam_dict = self.input()
        return bam_dict 

class extract_exon_annotation(luigi.Task):
    eisa_dir = '%s/%s/eisa' % (parameters().exp_dir, self.sample)
    
    def run(self):
        if not os.path.exists(self.eisa_dir):
            os.makedirs(self.eisa_dir)
        with open('%s/exons.gtf' % self.eisa_dir, 'a') as f:
            for line in open(parameters().genome_gtf):
                if "##" in line:
                    print(line, file = f)
                else:
                    cols = line.split()
                    if cols[2] == 'exon':
                        print(line, file = f)
        #awk_command = ["awk", "'$3 == "exon"'", "%s" % parameters().genome_gtf, ">", "%s/exons.gtf" % eisa_dir]
        #subprocess.call(awk_command)

    def output(self):
        return luigi.LocalTarget('%s/exons.gtf' % self.eisa_dir)

class filter_unassigned_mapped_reads(luigi.Task):
    
    sample = luigi.Parameter()
    bam_file_path = luigi.Parameter()
    exon_gtf = luigi.Parameters()
    eisa_dir = '%s/%s/eisa' % (parameters().exp_dir, 'test')

    def requires(self):
        return all_star_align(), extract_exon_annotation()

    def run(self):        
        intersect_command = [
                           'bedtools', 
                           'intersect', 
                           '-a %s' % self.bam_file_path,
                           '-b %s' % self.exon_gtf,
                           '-v', 
                           '>', 
                           '%s.intron.bam' % self.sample
                           ]
        subprocess.call(intersect_command)     
        
    def output(self):
        return luigi.LocalTarget('%s/%s.intron.bam' % (self.eisa_dir, self.sample))
   

class all_filter_unassigned_mapped_reads(luigi.Task):

    def requires(self):
        bam_dict = all_star_align().output() 
        return {
                s:filter_unassigned_mapped_reads(sample = s, bam_file_path = b.path) 
                    for s,b in bam_dict.items()
                }             

    def output(self):
        return self.input()

"""


class featureCounts(luigi.Task):
    
    sample = luigi.Parameter()
    bam_file = luigi.Parameter()    

    def requires(self):
        return all_star_align()

    def run(self):
        fC_dir = '%s/%s/featureCounts' % (parameters().exp_dir, self.sample)
        if not os.path.exists(fC_dir):
            os.makedirs(fC_dir)
        featureCounts_command = [
                                'featureCounts',
                                    '-F GTF',
                                    '-T %d' % cores,
                                    '-s %d' % stranded,
                                    '-a %s' % genome_gtf,
                                    '-o %s/%s.prelim.counts' % (fC_dir, self.sample),
                                    self.bam_file.path
                                ]
        for line in parameters().featureCounts.splitlines():
            featureCounts_command.append(line)
        fC = subprocess.Popen(featureCounts_command)
        fC.wait()
        if fC.returncode != 0:
            subprocess.call(['rm', '-rf', fC_dir])
        else:
            os.rename('%s/%s.prelim.counts' % (fC_dir, self.sample), 
                        '%s/%s.counts' % (fC_dir, self.sample))
            os.rename('%s/%s.prelim.counts.summary' % (fC_dir, self.sample), 
                        '%s/%s.counts.summary' % (fC_dir, self.sample))

    def output(self):
        fC_dir = '%s/%s/featureCounts' % (parameters().exp_dir, self.sample)
        return luigi.LocalTarget('%s/%s.counts' % (fC_dir, self.sample))


class all_featureCounts(luigi.Task):
    '''requires splits out samples and run featureCounts concurrently
    output returns dictionary with sample name and count file
    '''
    
    def requires(self):
        bam_dict = all_star_align().output() 
        return {
                s:featureCounts(sample = s, bam_file = b.path) 
                    for s,b in bam_dict.items()
                }             

    def output(self):
        counts_dict = self.input()
        return counts_dict

class luigi_count_matrix_postgres(luigi.postgres.CopyToTable):
    
    host = parameters().postgres_host
    database = parameters().postgres_database
    user = parameters().postgres_user
    password = parameters().postgres_password
    table = parameters().exp_name
        
    columns = [("Gene", "TEXT")]
    columns += [(name, "INT") for name in all_featureCounts().output()]                  

    def requires(self):
        return all_featureCounts()
    
    def rows(self): 
        count_files = [self.input()[y].path for y in self.input()]
        pandas_files = [
                        pd.read_table(self.input()[name].path, 
                            skiprows = 2,
                            index_col = 0,
                            names = ['Gene', 'Chr', 'Start', 'End', 
                                        'Strand', 'Length', name],
                            usecols = ['Gene', name],
                            header = None)
                        for name in self.input()
                        ]
        count_table = pd.concat(pandas_files, axis = 1).sort_index(axis=1)
        count_table = count_table.to_records()
        
        for row in count_table:
            yield(row)
"""

if __name__ == '__main__':
    luigi.run()
