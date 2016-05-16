#!/hpc/users/hagenj02/luigi_pipeline/vluigi/bin/python3.5
'''#!/home/jake/.virtualenvs/luigi/bin/python
    for testing on local machine
'''
import luigi
import luigi.postgres
#  import psycopg2
import subprocess
import os
import pandas as pd
import datetime
from sqlalchemy import create_engine
from sqlalchemy.schema import CreateSchema


class parameters(luigi.Config):
    '''
    Class to contain all parameters. Most of these will need to be set from
    luigi.cfg file in working directory, which is good because it acts as a log
    for parameters of experiment
    '''

    fastqs = luigi.Parameter(default=None)
    exp_dir = luigi.Parameter(default=os.getcwd())
    genome_fasta = luigi.Parameter(default=None)
    genome_gtf = luigi.Parameter(default=None)
    star_genome_folder = luigi.Parameter(default=None)
    read_length = luigi.IntParameter(default=100)
    stranded = luigi.BoolParameter(default=False)
    paried = luigi.IntParameter(default=0)
    cores = luigi.IntParameter(default=6)
    exp_name = luigi.Parameter(default=
                               datetime.date.today().strftime("%B%d,%Y"))


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


class fastqc(luigi.Task):
    '''Runs fastqc, quality control program on fastqs
    '''
    sample = luigi.Parameter()

    def require(self):
        return fastqs(self.sample)

    def output_dir(self):
        return '%s/fastqc' % parameters().exp_dir

    def run(self):
        try:
            os.makedirs('%s/fastqc' % self.output_dir())
        except OSError:
            pass

        fastqc_command = ['fastqc', self.input().path,
                          '-o', self.output_dir()]
        subprocess.call(fastqc_command)

    def output(self):
        file_name = os.path.basename(self.input().path).split(".")[0]
        return luigi.LocalTarget('%s/%s_fastqc.html' %
                                (self.output_dir(), file_name))


class star_index(luigi.Task):
    '''Index genome to be used with STAR aligner
    - More options can be passed to STAR by adding paramters to the config file
    '''
    star_index = luigi.Parameter(default="", significant=False)

    def run(self):
        try:
            os.makedirs(parameters().star_genome_folder)
        except OSError:
            pass

        star_command = ['STAR',
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

       More options can be passed to STAR by adding parameters to the
       config file
       The config file already contains default aligner options
    '''
    sample = luigi.Parameter()
    star_align = luigi.Parameter(default="", significant=False)

    def requires(self):
        return star_index(), fastqc(sample=self.sample)

    def output_dir(self):
        return '%s/%s/star' % (parameters().exp_dir, self.sample)

    def run(self):
        try:
            os.makedirs(self.output_dir())
        except OSError:
            pass

        star_command = ['STAR',
                        '--genomeDir %s' % parameters().star_genome_folder,
                        '--readFilesIn %s' % self.input()[1].path,
                        '--runThreadN %d' % parameters().cores,
                        '--outFileNamePrefix %s/%s.' %
                        (self.output_dir(), self.sample)
                        ]
        for line in self.star_align.splitlines():
            star_command.append(line)
        subprocess.call(star_command)

        if not os.path.isfile('%s/%s.Aligned.sortedByCoord.out.bam' %
                             (self.output_dir(), self.sample)):
            raise OSError("STAR could not create %s bam file" % self.sample)

    def output(self):
        return luigi.LocalTarget('%s/%s.Aligned.sortedByCoord.out.bam' %
                                (self.output_dir(), self.sample))


class gene_counter(luigi.Task):
    '''Use featureCounts from the subread package to count alignments
       that overlap a gene.
       This also serves as a base class for counting on
       different features besides genes
    '''
    sample = luigi.Parameter()
    annotation = luigi.Parameter(default=parameters().genome_gtf)
    bam_generator = luigi.TaskParameter(default=star_align)
    feature_to_count = luigi.Parameter(default='exon')
    grouper = luigi.Parameter(default='gene_id')
    feature_level = luigi.Parameter(default="")
    output_name = luigi.Parameter(default="gene")

    def output_dir(self):
        return '%s/%s/counts' % (parameters().exp_dir, self.sample)

    def requires(self):
        try:
            return self.bam_generator(sample=self.sample), require()
        except NameError:
            return self.bam_generator(sample=self.sample)

    def run(self):
        try:
            os.makedirs(self.output_dir())
        except OSError:
            pass

        bam_file = self.bam_generator(sample=self.sample).output().path
        featureCounts_command = ['featureCounts',
                                 '-T', '%d' % parameters().cores,
                                 '-t', '%s' % self.feature_to_count,
                                 '-g', '%s' % self.grouper, self.feature_level,
                                 '-o', '%s/%s.%s.counts' %
                                 (self.output_dir(),
                                     self.sample,
                                     self.output_name),
                                 '-a', self.annotation,
                                 bam_file]
        subprocess.call(featureCounts_command)

    def output(self):
        return luigi.LocalTarget('%s/%s.%s.counts' % (self.output_dir(), self.sample, self.output_name))


class exon_counter(gene_counter):
    feature_to_count = 'exon'
    feauture_level = '-f'
    grouper = 'transcript_id'
    output_name = 'exon'


class extract_exon_annotation(luigi.Task):
    '''extract only exons from annotation file,
       will be used to filter reads that align to an intron
    '''
    eisa_dir = '%s/eisa' % parameters().exp_dir

    def run(self):
        try:
            os.makedirs(self.eisa_dir)
        except OSError:
            pass

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
       Should end up with a bam file with reads that align to an intron region
       or do not align at all
    '''
    sample = luigi.Parameter()
    exon_gtf = luigi.Parameter(extract_exon_annotation().output().path)

    def output_dir(self):
        return '%s/eisa/%s' % (parameters().exp_dir, self.sample)

    def requires(self):
        return star_align(self.sample), extract_exon_annotation()

    def run(self):
        try:
            os.makedirs(self.output_dir())
        except OSError:
            pass

        bam_file = self.input()[0]
        intersect_command = [
            'bedtools', 'intersect', '-a', bam_file.path, '-b', self.exon_gtf, '-v'
                                    ]
        with open('%s/%s.intron.bam' % (self.output_dir(), self.sample), 'w') as f:
            subprocess.call(intersect_command, stdout=f)

    def output(self):
        return luigi.LocalTarget('%s/%s.intron.bam' % (self.output_dir(), self.sample))


class intron_counter(gene_counter):
    bam_generator = filter_nonexon
    feature_to_count = 'gene'
    output_name = 'intron'


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
    require = extract_protein_coding_annotation()
    output_name = "gene_protein_code"


class protein_coding_gene_intron_counter(gene_counter):
    require = extract_protein_coding_annotation
    bam_generator = luigi.TaskParameter(filter_nonexon)
    feature_to_count = 'gene'
    output_name = "intron_protein_code"
    annotation = extract_protein_coding_annotation().output().path


class multiQC_test(luigi.Task):

    def run(self):
        subprocess.call(['multiqc', parameters().exp_dir])


class all_counters(luigi.WrapperTask):
    def requires(self):
        print(fastqs(sample = "").run())
        yield {sample:gene_counter(sample = sample) for sample in fastqs(sample = "").run()}
        yield {sample:exon_counter(sample = sample) for sample in fastqs(sample = "").run()}
        yield {sample:intron_counter(sample = sample) for sample in fastqs(sample = "").run()}
        yield {sample:protein_coding_gene_counter(sample = sample) for sample in fastqs(sample = "").run()}
        yield {sample:protein_coding_gene_intron_counter(sample = sample) for sample in fastqs(sample = "").run()}


class postgres_count_matrix(luigi.Task):
    password = luigi.Parameter(significant=False)
    host = luigi.Parameter(significant=False)
    database = 'rna'
    user = luigi.Parameter(default='rna', significant=False)
    table = luigi.Parameter(default='gene_counts')
    feature_counter = luigi.TaskParameter(default=gene_counter, significant=False)

    def requires(self):
        return {x:self.feature_counter(sample=x) for x in fastqs(sample='').run()}

    def run(self):
        engine = create_engine('postgresql://%s:%s@%s/%s' %
                               (self.user, self.password, self.host, self.database))

        try:
            engine.execute(CreateSchema(parameters().exp_name))
        except:    # should catch psycopg2.ProgrammingError, but doesnt work
            pass

        pandas_files = [
                        pd.read_table(self.input()[name].path,
                                      skiprows=2,
                                      index_col=0,
                                      names=['Gene', 'Chr', 'Start', 'End',
                                             'Strand', 'Length', name],
                                      usecols=['Gene', name],
                                      header=None)
                        for name in self.input()
                        ]
        count_table = pd.concat(pandas_files, axis=1).sort_index(axis=1)
        count_table.to_csv("%s/%s.csv" % (parameters().exp_dir, self.table))
        count_table.to_sql(self.table, con=engine, schema=parameters().exp_name)

        # Taken from luigi source code, makes marker table and adds entry
        self.output().create_marker_table()
        connection = self.output().connect()
        self.output().touch(connection)
        connection.commit()
        connection.close()

    def output(self):
        return luigi.postgres.PostgresTarget(host=self.host,
                                             database=self.database,
                                             user=self.user,
                                             password=self.password,
                                             table=self.table,
                                             update_id=parameters().exp_name +
                                             '_' + self.table)


class all_count_matrix(luigi.WrapperTask):
    password = luigi.Parameter(significant=False)
    host = luigi.Parameter(significant=False)

    def requires(self):
        yield postgres_count_matrix(password=self.password, host=self.host)
        yield postgres_count_matrix(table="exon_counts",
                                    feature_counter=exon_counter,
                                    password=self.password, host=self.host)
        yield postgres_count_matrix(table="protein_gene_counts",
                                    feature_counter=protein_coding_gene_counter,
                                    password=self.password, host=self.host)
        yield postgres_count_matrix(table="intron_counts",
                                    feature_counter=intron_counter,
                                    password=self.password, host=self.host)
        yield postgres_count_matrix(table="protein_intron_counts",
                                    feature_counter=protein_coding_gene_intron_counter,
                                    password=self.password, host=self.host)


if __name__ == '__main__':
    luigi.run()
