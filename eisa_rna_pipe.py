import luigi
import os
import subprocess
import basic_rna_pipe as brp


class extract_exon_annotation(luigi.Task):
    '''extract only exons from annotation file,
       will be used to filter reads that align to an intron
    '''
    eisa_dir = '%s/eisa' % brp.configs().exp_dir

    def run(self):
        try:
            os.makedirs(self.eisa_dir)
        except OSError:
            pass

        with open('%s/exons.gtf' % self.eisa_dir, 'a') as f:
            for line in open(brp.configs().genome_gtf):
                if "##" in line:
                    print(line, end='', file=f)
                else:
                    cols = line.split()
                    if cols[2] == 'exon':
                        print(line, end='', file=f)

    def output(self):
        return luigi.LocalTarget('%s/exons.gtf' % self.eisa_dir)


class filter_nonexon(luigi.Task):
    '''Make bam file with reads that do not overlap any exons
       Should end up with a bam file with reads that align to an intron region
       or do not align at all
    '''
    sample = luigi.Parameter()

    def output_dir(self):
        return '%s/eisa/%s' % (brp.configs().exp_dir, self.sample)

    def requires(self):
        return brp.star_align(self.sample), extract_exon_annotation()

    def run(self):
        try:
            os.makedirs(self.output_dir())
        except OSError:
            pass

        bam_file = self.input()[0]
        exon_gtf = self.input()[1]
        intersect_command = ['bedtools', 'intersect',
                             '-a', bam_file.path,
                             '-b', exon_gtf.path,
                             '-v']
        with open('%s/%s_intron.bam' % (self.output_dir(), self.sample), 'w') as f:
            subprocess.call(intersect_command, stdout=f)

    def output(self):
        return luigi.LocalTarget('%s/%s_intron.bam' %
                                (self.output_dir(), self.sample))


class intron_counter(brp.gene_counter):
    bam_generator = filter_nonexon
    feature_to_count = 'gene'
    output_name = 'intron'


class nascent_count_matrix(brp.postgres_count_matrix):
    table = "intron_counts"
    feature_counter = intron_counter
