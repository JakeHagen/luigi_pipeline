import luigi
import subprocess
import os

genome_fasta = '/sc/orga/projects/Houton_Sander/genomes/rna_star_99_gencode_human_8-17-15/GRCh38.primary_assembly.genome.fa'
genome_gtf = '/sc/orga/projects/Houton_Sander/genomes/rna_star_99_gencode_human_8-17-15/gencode.v23.annotation.gtf'
star_genome_folder = '/sc/orga/projects/Houton_Sander/genomes/rna_star_99_gencode_human_8-17-15/star_genome'
read_length = 100
cores = 4
genome_location = '/sc/orga/projects/Houton_Sander/genomes/rna_star_99_gencode_human_8-17-15/star_genome'
fastq_dictionary = {s01:'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi1_CAGATC_L008_R1_001.C6673ACXX.fastq.gz',
					s02:'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi2_ATCACG_L008_R1_001.C6673ACXX.fastq.gz',
					s03:'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi3_TCGGCA_L008_R1_001.C6673ACXX.fastq.gz',
					c13:'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lc13_CGATGT_L008_R1_001.C6673ACXX.fastq.gz',
					c14:'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lc14_GATCAG_L008_R1_001.C6673ACXX.fastq.gz',
					c15:'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lc15_CTTGTA_L008_R1_001.C6673ACXX.fastq.gz'  
					}

class genomeFiles(luigi.Task):
	def output(self):
		return luigi.LocalTarget(genome_fasta),  luigi.LocalTarget(genome_gtf)


class index_STAR_genome(luigi.Task):
	def requires(self):
		return genomeFiles()
	
	def run(self):
		if not os.path.exists(star_genome_folder):
			os.mkdir(star_genome_folder))
		star_command = [
						'STAR', '--runThreadN %d' % cores, '--runMode genomeGenerate', 
								'--genomeDir %s' % star_genome_folder, '--genomeFastaFiles %s' % genome_fasta, 
								'--sjdbGTFfile %s' % genome_gtf, 'sjdbOverhang %d' % (read_length - 1)
						]
		subprocess.call(star_command)
	
	def output(self):
		return luigi.LocalTarget(star_genome_folder)


class fastqs(luigi.Task):
    def output(self):
        return fastq_dictionary


class star_align(luigi.Task):
	def requires(self):
		return index_STAR_genome(), fastqs()
	
	def run(self):
		for sample,path in fastq_dictionary:
			star_command = [
							'STAR', '--runThreadN %d' % cores, '--runMode
		

if __name__ == '__main__':
	luigi.run()


