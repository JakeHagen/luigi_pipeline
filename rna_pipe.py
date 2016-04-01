#!/home/jake/.local/share/virtualenvs/luigi/bin/python3.5
#####!/hpc/users/hagenj02/luigi_pipeline/vluigi/bin/python3.5

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
    '''Takes fastqs from parameters (specified in python.cfg),
       makes dictionary, and returns fastq based on sample
    '''
    sample = luigi.Parameter()

    def run(self):
        fastq_dict = {}
        for line in parameters().fastqs.splitlines():
            sample, path = line.split(":")
            fastq_dict[sample] = path
        return fastq_dict

    def output(self):
        fastq_dict = self.run()
        fastq_file_path = fastq_dict[self.sample]
        return luigi.LocalTarget(fastq_file_path)


class star_index(luigi.Task):
    '''Index genome to be used with STAR aligner

       More options can be passed to STAR by adding paramters to the config file
    '''
    star_index = luigi.Parameter(default = "", significant = False)

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
        if not os.path.isfile('%s/Genome' % parameters().star_genome_folder):
            raise OSError("star_index/Genome file could not be created")

    def output(self):
        return luigi.LocalTarget('%s/Genome' % parameters().star_genome_folder)


class star_align(luigi.Task):
    '''Align fastq to previously indexed genome

       More options can be passed to STAR by adding parameters to the config file
       The config file already contains default aligner options
    '''
    sample = luigi.Parameter()
    star_align = luigi.Parameter(default = "", significant = False)

    def requires(self):
        return star_index(), fastqs(sample = self.sample)

    def run(self):
        if not os.path.exists('%s/%s/star' % (parameters().exp_dir, self.sample)):
            os.makedirs('%s/%s/star' % (parameters().exp_dir, self.sample))
        star_command = [
                        'STAR',
                        '--genomeDir %s' % parameters().star_genome_folder,
                        '--readFilesIn %s' % self.input()[1].path,
                        '--runThreadN %d' % parameters().cores,
                        '--outFileNamePrefix %s/%s/star/%s.' %
                            (parameters().exp_dir, self.sample, self.sample)
                        ]
        #Append the extra parameters we want from the config file
        for line in self.star_align.splitlines():
            star_command.append(line)
        subprocess.call(star_command)
        if not os.path.isfile('%s/%s/star/%s.Aligned.sortedByCoord.out.bam' % (parameters().exp_dir, self.sample, self.sample)):
            raise OSError("STAR bam could not be created")

    def output(self):
        return luigi.LocalTarget('%s/%s/star/%s.Aligned.sortedByCoord.out.bam' % (parameters().exp_dir, self.sample, self.sample))


class gene_counter(luigi.Task):
    '''Use featureCounts from the subread package to count alignments that overlap
       a gene.
       This also serves as a base class for counting on different features besides genes
    '''
    sample = luigi.Parameter()
    annotation = luigi.Parameter(default = parameters().genome_gtf)
    feature = luigi.Parameter(default = 'gene')
    require = luigi.TaskParameter(default = star_align)
    output_name = luigi.Parameter(default = 'gene')
    bam_generator = luigi.TaskParameter(default = star_align)

    def output_dir(self):
        return  '%s/%s/counts' % (parameters().exp_dir, self.sample)

    def requires(self):
        return self.require(sample = self.sample)

    def run(self):
        if not os.path.exists(self.output_dir()):
            os.makedirs(self.output_dir())

        featureCounts_command = ['featureCounts', '-T', '%d' % parameters().cores,
                                    '-t', 'gene', '-g', 'gene_id',
                                    '-o', '%s/%s.%s.counts' % (self.output_dir(), self.sample, self.output_name),
                                    '-a', self.annotation,
                                    self.bam_generator(sample = self.sample).output().path]
        subprocess.call(featureCounts_command)

    def output(self):
        return luigi.LocalTarget('%s/%s.%s.counts' % (self.output_dir(), self.sample, self.output_name))


class exon_counter(gene_counter):
    feature = 'exon'
    output_name = 'exon'


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
    '''Make bam file with reads that do not overlap any exons
    To be used to count reads that map to introns
    '''
    sample = luigi.Parameter()
    exon_gtf = luigi.Parameter(extract_exon_annotation().output().path)

    def requires(self):
        return star_align(self.sample), extract_exon_annotation()

    def run(self):
        eisa_dir = '%s/eisa/%s' % (parameters().exp_dir, self.sample)
        bam_file = self.input()[0]

        if not os.path.exists(eisa_dir):
            os.makedirs(eisa_dir)

        intersect_command = [
            'bedtools', 'intersect', '-a', '%s' % bam_file.path, '-b', '%s' % self.exon_gtf, '-v'
                           ]

        with open('%s/%s.intron.bam' % (eisa_dir, self.sample), 'w') as f:
            subprocess.call(intersect_command, stdout=f)

    def output(self):
        eisa_dir = '%s/eisa/%s' % (parameters().exp_dir, self.sample)
        return luigi.LocalTarget('%s/%s.intron.bam' % (eisa_dir, self.sample))


class intron_counter(gene_counter):
    require = luigi.TaskParameter(filter_nonexon)
    feature = 'intron'
    output_name = 'intron'
    bam_generator = luigi.TaskParameter(default = filter_nonexon)


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
    annotation = extract_protein_coding_annotation().output().path
    output_name = 'protein_coding'

    def requires(self):
        return star_align(self.sample), extract_protein_coding_annotation()

class protein_coding_gene_intron_counter(gene_counter):
    annotation = extract_protein_coding_annotation().output().path
    output_name = 'protein_coding_intron'

    def requires(self):
        return filter_nonexon(self.sample), extract_protein_coding_annotation()

def table_maker(feature):
    return parameters().exp_name + '_' + feature

class postgres_gene_count_matrix(luigi.postgres.CopyToTable):
    password = luigi.Parameter(significant = False)
    host = luigi.Parameter(significant = False)
    database = luigi.Parameter(default = 'rna', significant = False)
    user = luigi.Parameter(default = 'rna', significant = False)
    #feature = luigi.Parameter(default = "gene")
    table = luigi.Parameter(default = parameters().exp_name + '_' + 'gene_counts')
    feature_counter = luigi.TaskParameter(default = gene_counter, significant = False)

    columns = [("Gene", "TEXT")]
    columns += [(name, "INT") for name in fastqs(sample = '').run()]

    def requires(self):
        return {x:self.feature_counter(sample = x) for x in fastqs(sample = '').run()}

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

class postgres_exon_count_matrix(postgres_gene_count_matrix):
    table = luigi.Parameter(default = parameters().exp_name + '_' + 'exon_counts')
    feature_counter = exon_counter

class postgres_protein_genes_count_matrix(postgres_gene_count_matrix):
    table = luigi.Parameter(default = parameters().exp_name + '_' + 'protein_genes_counts')
    feature_counter = protein_coding_gene_counter

class postgres_intron_count_matrix(postgres_gene_count_matrix):
    table = luigi.Parameter(default = parameters().exp_name + '_' + 'intron_counts')
    feature_counter = intron_counter

class postgres_protein_intron_count_matrix(postgres_gene_count_matrix):
    table = luigi.Parameter(default = parameters().exp_name + '_' + 'intron_protein_counts')
    feature_counter = protein_coding_gene_intron_counter

class all_count_matrix(luigi.WrapperTask):
    password = luigi.Parameter(significant = False)
    def requires(self):
        yield postgres_gene_count_matrix(self.password)
        yield postgres_exon_count_matrix(self.password)
        yield postgres_protein_genes_count_matrix(self.password)
        yield postgres_intron_count_matrix(self.password)
        yield postgres_protein_intron_count_matrix(self.password)

if __name__ == '__main__':
    luigi.run()
