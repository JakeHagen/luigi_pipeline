#!/hpc/users/hagenj02/luigi_pipeline/vluigi/bin/python3.5

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


class star_index(luigi.Task):
    
    star_index = luigi.Parameter(default = "")
        
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
        for line in self.star_index.splitlines():
            star_command.append(line)
        subprocess.call(star_command)
            
    def output(self):
        return luigi.LocalTarget('%s/Genome' % parameters().star_genome_folder)

class split_group(luigi.Task):
    input_task = luigi.TaskParameter()
    task = luigi.TaskParameter()

    def requires(self):
        return{s:self.task(sample = s, luigi_file = f) for s,f in self.input_task().output().items()}

    def output(self):
        return self.input()
        


class star_align(luigi.Task):
    
    sample = luigi.Parameter()
    luigi_file = luigi.Parameter()
    star_align = luigi.Parameter(default = "")
    
    def requires(self):
        return star_index(), fastqs()

    def run(self):
        if not os.path.exists('%s/%s/star' % (parameters().exp_dir, self.sample)):
            os.makedirs('%s/%s/star' % (parameters().exp_dir, self.sample))
        star_command = [
                        'STAR',
                        '--genomeDir %s' % parameters().star_genome_folder,
                        '--readFilesIn %s' % self.luigi_file.path,
                        '--runThreadN %d' % parameters().cores,
                        '--outFileNamePrefix %s/%s/star/%s.' %
                            (parameters().exp_dir, self.sample, self.sample)
                        ]
        #Append the extra parameters we want from the config file
        for line in self.star_align.splitlines():
            star_command.append(line)
        subprocess.call(star_command)
        os.rename('%s/%s/star/%s.Aligned.sortedByCoord.out.bam' % (parameters().exp_dir, self.sample, self.sample), 
                    '%s/%s/star/%s.bam' % (parameters().exp_dir, self.sample, self.sample))
    def output(self):
        return luigi.LocalTarget('%s/%s/star/%s.bam' 
                                    % (parameters().exp_dir, self.sample, self.sample))


class sg_star_align(split_group):
    input_task = luigi.TaskParameter(fastqs)
    task = luigi.TaskParameter(star_align)

 
class gene_counter(luigi.Task):

    sample = luigi.Parameter()
    luigi_file = luigi.Parameter()
    annotation = luigi.Parameter(default = parameters().genome_gtf)
    feature = luigi.Parameter(default = 'gene')
    require = luigi.TaskParameter(sg_star_align)
    
    def output_dir(self):
        return  '%s/%s/counts' % (parameters().exp_dir, self.sample)

    def requires(self):
        return self.require()

    def run(self):
        if not os.path.exists(self.output_dir()):
            os.makedirs(self.output_dir())
        featureCounts_command = ['featureCounts', '-T', '%d' % parameters().cores,
                                    '-t', self.feature, '-g', 'gene_id',
                                    '-o', '%s/%s.%s.counts' % (self.output_dir(), self.sample, self.feature),
                                    '-a', self.annotation,
                                    self.luigi_file.path]

        subprocess.call(featureCounts_command)

    def output(self):
        return luigi.LocalTarget('%s/%s.%s.counts' % (self.output_dir(), self.sample, self.feature))

class sg_gene_counter(split_group):
    input_task = luigi.TaskParameter(sg_star_align)
    task = luigi.TaskParameter(gene_counter)

class exon_counter(luigi.Task):
    feature = 'exon'

class sg_exon_counter(split_group):
    input_task = luigi.TaskParameter(sg_star_align)
    task = luigi.TaskParameter(exon_counter)


class extract_exon_annotation(luigi.Task):

    eisa_dir = '%s/eisa' % parameters().exp_dir

    def run(self):
        if not os.path.exists(self.eisa_dir):
            os.makedirs(self.eisa_dir)
        
        with open('%s/exons.gtf' % self.eisa_dir, 'a') as f:
            for line in open(parameters().genome_gtf):
                if "##" in line:
                    print(line, end = "", file = f)
                else:
                    cols = line.split()
                    if cols[2] == 'exon':
                        print(line, end = "", file = f)

    def output(self):
        return luigi.LocalTarget('%s/exons.gtf' % self.eisa_dir)

class filter_nonexon(luigi.Task):

    sample = luigi.Parameter()
    bam_file_path = luigi.Parameter()
    exon_gtf = extract_exon_annotation().output().path

    def requires(self):
        return split_group_star_align(), extract_exon_annotation()

    def run(self):
        eisa_sample_dir = '%s/eisa/%s' % (parameters().exp_dir, self.sample)
        if not os.path.exists(eisa_sample_dir):
            os.makedirs(eisa_sample_dir)
        
        intersect_command = [
            'bedtools', 'intersect', '-a', '%s' % self.bam_file_path, '-b', '%s' % self.exon_gtf, '-v'
                           ]
        with open('%s/%s.intron.bam' % (eisa_sample_dir, self.sample), 'w') as f:
            subprocess.call(intersect_command, stdout=f)

    def output(self):
        eisa_sample_dir = '%s/eisa/%s' % (parameters().exp_dir, self.sample)
        return luigi.LocalTarget('%s/%s.intron.bam' % (eisa_sample_dir, self.sample))


class sg_filter_nonexon(split_group):
    input_task = luigi.TaskParameter(sg_star_align)
    task = luigi.TaskParameter(filter_nonexon)


class intron_counter(gene_counter):
    require = luigi.TaskParameter(sg_filter_nonexon)


class sg_intron_counter(split_group):
    input_task = luigi.TaskParameter(sg_filter_nonexon)
    task = luigi.TaskParameter(intron_counter)


class extract_protein_coding_annotation(luigi.Task):

    def run(self):
        with open('%s/protein_coding.gtf' % parameters().exp_dir, 'a') as f:
            for line in open(parameters().genome_gtf):
                if "##" in line:
                    print(line, end = "", file = f)
                else:
                    if 'protein_coding' in line:
                        print(line, end = "", file = f)

    def output(self):
        return luigi.LocalTarget('%s/protein_coding.gtf' % parameters().exp_dir)


class protein_coding_gene_counter(gene_counter):
    require = luigi.TaskParameter(extract_protein_coding_annotation)
    annotation =luigi.Parameter(extract_protein_coding_annotation().output().path)
    

class sg_protein_coding_gene_counter(split_group):
    input_task = luigi.TaskParameter(sg_star_align)
    task = luigi.TaskParameter(protein_coding_gene_counter)


class sg_protein_coding_gene_intron_counter(split_group):
    input_task = luigi.TaskParameter(sg_filter_nonexon)
    task = luigi.TaskParameter(protein_coding_gene_counter)


class postgres_count_matrix(luigi.Task):#postgres.CopyToTable):

    #password = luigi.Parameter(significant = False)
    #host = luigi.Parameter(default = 'localhost')
    #database = luigi.Parameter(default = 'RNA')
    #user = luigi.Parameter(default = 'hagenj02')
    #feature = luigi.Parameter(default = "gene")
    #table = luigi.Parameter(default = parameters().exp_name)
    #feature_counter = luigi.TaskParameter(default = gene_counter)
    
    def cols(self):
        columns = [("Gene", "TEXT")]
        columns += [(name, "INT") for name in gene_counter().output()]
        return columns
    
    columns = self.cols()

    def run(self):
        print(self.columns)
"""    
    def requires(self):
        return self.feature_counter()

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
        count_table.to_csv("%s/%s.csv" % (parameters().exp_dir, self.table))
        count_table = count_table.to_records()

        for row in count_table:
            yield(row)
"""

if __name__ == '__main__':
    luigi.run()



