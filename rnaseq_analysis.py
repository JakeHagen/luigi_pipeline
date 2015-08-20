import luigi
import subprocess
import os

#will output alignment and count data to the 'working_dir' but the genome index will be created and stored in 'star_genome_folder'
working_dir = '/sc/orga/projects/argmac01a'
genome_fasta = '/sc/orga/projects/Houton_Sander/genomes/rna_star_99_gencode_human_8-17-15/GRCh38.primary_assembly.genome.fa'
genome_gtf = '/sc/orga/projects/Houton_Sander/genomes/rna_star_99_gencode_human_8-17-15/gencode.v23.annotation.gtf'
star_genome_folder = '/sc/orga/projects/Houton_Sander/genomes/rna_star_99_gencode_human_8-17-15/star_genome'
read_length = 100
cores = 12 
fastq_dictionary = {
					's01':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi1_CAGATC_L008_R1_001.C6673ACXX.fastq.gz',
					's02':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi2_ATCACG_L008_R1_001.C6673ACXX.fastq.gz',
					's03':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi3_TCGGCA_L008_R1_001.C6673ACXX.fastq.gz',
					'c13':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lc13_CGATGT_L008_R1_001.C6673ACXX.fastq.gz',
					'c14':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lc14_GATCAG_L008_R1_001.C6673ACXX.fastq.gz',
					'c15':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lc15_CTTGTA_L008_R1_001.C6673ACXX.fastq.gz'  
					}

class genomeFiles(luigi.Task):
	def output(self):
		return luigi.LocalTarget(genome_fasta),  luigi.LocalTarget(genome_gtf)

class index_STAR_genome(luigi.Task):
	def requires(self):
		return genomeFiles()
	
	def run(self):
		if not os.path.exists(star_genome_folder):
			os.mkdir(star_genome_folder)
		subprocess.call(['cd %s/..' % star_genome_folder])
		star_command = [
						'STAR', 
								'--runThreadN %d' % cores, 
								'--runMode genomeGenerate', 
								'--genomeDir %s' % star_genome_folder, 
								'--genomeFastaFiles %s' % genome_fasta, 
								'--sjdbGTFfile %s' % genome_gtf, 
								'--sjdbOverhang %d' % (read_length - 1)
						]
		star = subprocess.Popen(star_command)
		star.wait()
		if star.returncode == 0:
			open('%s.token' % star_genome_folder, 'a').close() 
	def output(self):
		return luigi.LocalTarget('%s.token' % star_genome_folder)



class make_star_align_dir(luigi.Task):
	def requires(self):
		return index_STAR_genome() 
	def run(self):
		os.mkdir('%s/star_align' % working_dir)
	def output(self):
		return luigi.LocalTarget('%s/star_align' % working_dir)

class all_star_align(luigi.WrapperTask):
	def requires(self):
		for s,p in fastq_dictionary.items():
			yield star_align(sample = s, path = p)

class star_align(luigi.Task):
	sample = luigi.Parameter()
	path = luigi.Parameter()
	
	def requires(self):
		return index_STAR_genome(), make_star_align_dir()
	
	def run(self):
		if not os.path.exists('%s/star_align/%s' % (working_dir,self.sample)):
			os.mkdir('%s/star_align/%s' % (working_dir,self.sample))	
		
		s_command = [
					'STAR', 
						'--genomeDir %s' % star_genome_folder, 
						#'--sjdbGTFfile %s' % genome_gtf, 
						'--readFilesIn %s' % self.path, 
						'--readFilesCommand zcat', 
						'--runThreadN %d' % cores, #(cores/len(fastq_dictionary)), 
						'--outSAMmode Full', 
						'--outReadsUnmapped Fastx', 	
						'--chimSegmentMin 15', 
						'--chimJunctionOverhangMin 15', 	
						'--outSAMstrandField intronMotif', 
						'--outFilterType BySJout', 
						'--outFilterIntronMotifs RemovenonconicalUnannotated', 
						#'--genomeLoad LoadAndRemove',
						#'--limitBAMsortRAM 10000000000', 
						'--outSAMtype BAM SortedByCoordinate',
						'--outFileNamePrefix %s/star_align/%s/%s.' % (working_dir, self.sample, self.sample) 
					]
		star = subprocess.Popen(s_command)
		star.wait()
		if star.returncode == 0:
			open('%s/star_align/%s/star_align.token' % (working_dir, self.sample), 'a').close
	def output(self):
		return luigi.LocalTarget('%s/star_align/%s/star_align.token' % (working_dir, self.sample))
									
		

if __name__ == '__main__':
	luigi.run()

