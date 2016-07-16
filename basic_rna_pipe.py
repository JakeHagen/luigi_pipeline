import luigi
import luigi.postgres
import os
import subprocess
from sqlalchemy import create_engine
from sqlalchemy.schema import CreateSchema
import psycopg2
import pandas as pd


class configs(luigi.Config):
    '''
    Class to contain common parameters. Most of these can be set from
    luigi.cfg file in working directory, which is good because it acts as a log
    for parameters of experiment
    '''
    star_genome_folder = luigi.Parameter()
    exp_dir = luigi.Parameter(default=os.getcwd())
    genome_fasta = luigi.Parameter()
    genome_gtf = luigi.Parameter()
    read_length = luigi.IntParameter(default=100)
    stranded = luigi.BoolParameter(default=False)
    paried = luigi.IntParameter(default=False)
    cores = luigi.IntParameter(default=1)
    exp_name = luigi.Parameter()


class star_index(luigi.Task):
    '''Index genome to be used with STAR aligner

       More options can be passed to STAR by adding paramters to config file
    '''
    star_index_extra_params = luigi.Parameter(default='', significant=False)

    def run(self):
        try:
            os.makedirs(configs().star_genome_folder)
        except OSError:
            pass

        star_command = ['STAR',
                        '--runThreadN %d' % configs().cores,
                        '--runMode genomeGenerate',
                        '--genomeDir %s' % configs().star_genome_folder,
                        '--genomeFastaFiles %s' % configs().genome_fasta,
                        '--sjdbGTFfile %s' % configs().genome_gtf,
                        '--sjdbOverhang %d' % configs().read_length - 1,
                        ]
        for extra_param in self.star_index_extra_params.splitlines():
            star_command.append(extra_param)
        subprocess.call(star_command)

        if not os.path.isfile('%s/Genome' % configs().star_genome_folder):
            raise OSError("star_index/Genome file could not be created")

    def output(self):
        return luigi.LocalTarget('%s/Genome' % configs().star_genome_folder)


class fastqs(luigi.Task):
    '''Takes fastqs from parameters (specified in python.cfg),
       makes dictionary, and returns fastq based on sample

       fastqs.run() can also be used to generate dict with samples and path,
       useful for running all samples through one step in pipeline
    '''
    sample = luigi.Parameter(default=None)
    sample_fastqs = luigi.Parameter(significant=False)

    def output(self):
        fastq_dict = {}
        for line in self.sample_fastqs.splitlines():
            sample, path = line.split(":")
            try:
                path_one, path_two = path.split()
                path = path_one, path_two
            except ValueError:
                pass
            fastq_dict[sample] = path
        if self.sample is None:
            return fastq_dict
        elif type(fastq_dict[self.sample]) is tuple:
            return {'fastq_pair_1':
                    luigi.LocalTarget(fastq_dict[self.sample][0]),
                    'fastq_pair_2':
                    luigi.LocalTarget(fastq_dict[self.sample][1])}
        else:
            return luigi.LocalTarget(fastq_dict[self.sample])


class star_align(luigi.Task):
    '''Align fastq to previously indexed genome

       More options can be passed to STAR by adding parameters to the
       config file
       The config file already contains default aligner options
    '''
    sample = luigi.Parameter()
    star_align_extra_params = luigi.Parameter(default="", significant=False)

    def requires(self):
        return star_index(), fastqs(sample=self.sample)

    def output_dir(self):
        return '%s/star/%s' % (configs().exp_dir, self.sample)

    def run(self):
        try:
            os.makedirs(self.output_dir())
        except OSError:
            pass

        if type(self.input()[1]) is dict:
            fastqs = ''
            for x in self.input()[1]:
                fastqs += self.input()[1][x].path
                fastqs += ' '
        else:
            fastqs = self.input()[1].path
        star_command = ['STAR',
                        '--genomeDir %s' % configs().star_genome_folder,
                        '--readFilesIn %s' % fastqs,
                        '--runThreadN %d' % configs().cores,
                        '--outFileNamePrefix %s/%s.' %
                        (self.output_dir(), self.sample)
                        ]
        for line in self.star_align_extra_params.splitlines():
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
    annotation = luigi.Parameter(default=configs().genome_gtf)
    bam_generator = luigi.TaskParameter(default=star_align)
    feature_to_count = luigi.Parameter(default='exon')
    grouper = luigi.Parameter(default='gene_id')
    feature_level = luigi.Parameter(default='')
    output_name = luigi.Parameter(default='gene')

    def output_dir(self):
        return '%s/counts/%s' % (configs().exp_dir, self.sample)

    def requires(self):
        try:
            return self.bam_generator(sample=self.sample), self.require()
        except AttributeError:
            return self.bam_generator(sample=self.sample)

    # featureCounts command into function so it can be easily subclassed
    # for very specialized use cases
    def featureCounts_command(self):
        bam_file = self.bam_generator(sample=self.sample).output().path
        featCounts_command = ['featureCounts',
                              '-T', '%d' % configs().cores,
                              '-t', '%s' % self.feature_to_count,
                              '-g', '%s' % self.grouper, self.feature_level,
                              '-o', '%s/%s_%s.counts' %
                              (self.output_dir(),
                               self.sample,
                               self.output_name),
                              '-a', self.annotation,
                              bam_file]
        return featCounts_command

    def run(self):
        try:
            os.makedirs(self.output_dir())
        except OSError:
            pass
        subprocess.call(self.featureCounts_command())

    def output(self):
        return luigi.LocalTarget(
            '%s/%s_%s.counts' %
            (self.output_dir(), self.sample, self.output_name)
            )


class postgres_count_matrix(luigi.Task):
    '''insert compliled counts into postgres database
       Many generic parameters because this class will be used again
       to isert different count matrices into postgres
    '''
    password = luigi.Parameter(significant=False)
    host = luigi.Parameter(significant=False)
    database = luigi.Parameter(default='rna')
    user = luigi.Parameter(default='rna', significant=False)
    table = luigi.Parameter(default='gene_counts')
    feature_counter = luigi.TaskParameter(default=gene_counter,
                                          significant=False)

    def requires(self):
        return {samp: self.feature_counter(sample=samp)
                for samp in fastqs().output()}

    def run(self):
        engine = create_engine('postgresql://%s:%s@%s/%s' %
                               (self.user, self.password,
                                self.host, self.database)
                               )

        try:
            engine.execute(CreateSchema(configs().exp_name))
        except:    # should catch psycopg2.ProgrammingError, but doesnt work
            pass

        # compile counts files into one matrix
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

        # write matrix to csv and insert in postgres
        count_table.to_csv("%s/%s.csv" % (configs().exp_dir, self.table))
        count_table.to_sql(self.table, con=engine, schema=configs().exp_name)

        # Taken from luigi source code, makes marker table and adds entry
        # This is what lets luigi know task is already completed
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
                                             update_id=configs().exp_name +
                                             '_' + self.table)
