import luigi
import basic_rna_pipe as brp
import eisa_rna_pipe as erp


class exon_counter(brp.gene_counter):
    feature_to_count = 'exon'
    feauture_level = '-f'
    grouper = 'transcript_id'
    output_name = 'exon'


class extract_protein_coding_annotation(luigi.Task):

    def run(self):
        with open('%s/protein_coding.gtf' % brp.configs().exp_dir, 'a') as f:
            for line in open(brp.configs().genome_gtf):
                if "##" in line:
                    print(line, end='', file=f)
                else:
                    if 'protein_coding' in line:
                        print(line, end='', file=f)

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


# These classes could have been kept as is and given parameters instead of
# subclassing. I thought this was easier, with less parameters to specify
class protein_gene_count_matrix(brp.postgres_count_matrix):
    table = "protein_gene_counts"
    feature_counter = protein_coding_gene_counter


class protein_intron_count_matrix(brp.postgres_count_matrix):
    table = "protein_intron_counts"
    feature_counter = protein_coding_gene_intron_counter
