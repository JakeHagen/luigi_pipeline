import luigi
import os
import basic_rna_pipe as brp
import eisa_rna_pipe as erp


class exon_counter(brp.gene_counter):
    feature_to_count = 'exon'
    feauture_level = '-f'
    grouper = 'transcript_id'
    output_name = 'exon'


class extract_protein_coding_annotation(luigi.Task):
    
    def make_protein_coding_gtf(self):
        with open('%s/protein_coding.gtf' % brp.configs().exp_dir, 'a') as f:
            for line in open(brp.configs().genome_gtf):
                if "##" in line:
                    print(line, end='', file=f)
                else:
                    if 'protein_coding' in line:
                        print(line, end='', file=f)

    def run(self):
        try:
            self.make_protein_coding_gtf()
        except FileNotFoundError:
            os.makedirs(brp.configs().exp_dir)
            self.make_protein_coding_gtf()

    def output(self):
        return luigi.LocalTarget('%s/protein_coding.gtf' %
                                 brp.configs().exp_dir)


class protein_coding_gene_counter(brp.gene_counter):
    annotation = extract_protein_coding_annotation().output().path
    require = extract_protein_coding_annotation
    output_name = "gene_protein_code"


class protein_coding_gene_intron_counter(brp.gene_counter):
    require = extract_protein_coding_annotation
    bam_generator = luigi.TaskParameter(erp.filter_nonexon)
    feature_to_count = 'gene'
    output_name = "intron_protein_code"
    annotation = extract_protein_coding_annotation().output().path


class viral_counter(brp.gene_counter):

    def featureCounts_command(self):
        bam_file = self.bam_generator(sample=self.sample).output().path
        featCounts_command = ['featureCounts',
                              '-T', '%d' % configs().cores,
                              '-t', '%s' % self.feature_to_count,
                              '-g', '%s' % self.grouper, self.feature_level,
                              '-M',  # count multimappers
                              '-o', '%s/%s_%s.counts' %
                              (self.output_dir(),
                               self.sample,
                               self.output_name),
                              '-a', self.annotation,
                              bam_file]
        return featCounts_command

# These classes could have been kept as is and given parameters instead of
# subclassing. I thought this was easier, with less parameters to specify
class protein_gene_count_matrix(brp.postgres_count_matrix):
    table = "protein_gene_counts"
    feature_counter = protein_coding_gene_counter


class protein_intron_count_matrix(brp.postgres_count_matrix):
    table = "protein_intron_counts"
    feature_counter = protein_coding_gene_intron_counter
