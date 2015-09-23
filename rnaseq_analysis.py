import luigi
import subprocess
import os


#will output will be, can be pointed at an existing folder or the pipeline will create it
wkdir = luigi.Parameter(default = '/sc/orga/projects/argmac01a/luigi_rna_ensembl37')
# Only needed if STAR genome needs to be generated
genome_fasta = luigi.Parameter(default = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa') 
genome_gtf = luigi.Parameter(default = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/Homo_sapiens.GRCh37.75.gtf')

"""Either where the STAR genome is located or were it will be created
   IMPORTANT: if the STAR genome needs to be created, do not create the folder,
   Let the pipeline create it."""
star_genome_folder = luigi.Parameter(default = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/star_genome_ensembl37_99')

# RNAseq experiment info
read_length = luigi.FloatParameter(default = 100)
stranded = luigi.FloatParameter(default = 0)
#paried = ?
#prep = polyA or ribozero

# Should be multiple of samples
cores = luigi.FloatParameter(default = 6)

# Needs to be sample name (can be anything) and location of fastq.gz
fastq_dictionary = {'s01':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi1_CAGATC_L008_R1_001.C6673ACXX.fastq.gz',
                    's02':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi2_ATCACG_L008_R1_001.C6673ACXX.fastq.gz',
                    's03':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi3_TCGGCA_L008_R1_001.C6673ACXX.fastq.gz',
                    'c13':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lc13_CGATGT_L008_R1_001.C6673ACXX.fastq.gz',
                    'c14':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lc14_GATCAG_L008_R1_001.C6673ACXX.fastq.gz',
                    'c15':'/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lc15_CTTGTA_L008_R1_001.C6673ACXX.fastq.gz',
                    }
# Not really need, was originally created for testing
class genomeFiles(luigi.Task):
    def output(self):
        return (luigi.LocalTarget(genome_fasta), luigi.LocalTarget(genome_gtf))

class index_STAR_genome(luigi.Task):
    def requires(self):
        return genomeFiles()

    def run(self):
        if not os.path.exists(star_genome_folder):
            os.makedirs(star_genome_folder) 
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
        if star.returncode != 0:
            subprocess.call(['rm -rf %s' % star_genome_folder])
        os.rename('Log.out', '%s/Log.out' % star_genome_folder)
    def output(self):
        return luigi.LocalTarget('%s/Log.out' % star_genome_folder)



class wk_dir(luigi.Task):
    def requires(self):
        return index_STAR_genome()
    def run(self):
        os.makedirs(wkdir)
    def output(self):
        return luigi.LocalTarget(wkdir)

class star_align(luigi.Task):
    sample = luigi.Parameter()
    path = luigi.Parameter()
    def requires(self):
        return index_STAR_genome(), wk_dir()

    def run(self):
        os.makedirs('%s/%s/star' % (wkdir, self.sample))
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
                    '--outFilterIntronMotifs RemoveNoncanonicalUnannotated',
                    #'--genomeLoad LoadAndRemove',
                    #'--limitBAMsortRAM 15000000000',
                    '--outSAMtype BAM SortedByCoordinate',#Unsorted
                    '--outFileNamePrefix %s/%s/star/%s.' % (wkdir, self.sample, self.sample)
                    ]
        star = subprocess.Popen(s_command)
        star.wait()
        # Below removes folder if star command was unsuccessful 
        # (files remain when star fails)
        # renaming bam insures that if batch job gets cut while star is running
        # the task will not be marked as complete 
        if star.returncode != 0:
            subprocess.call(['rm', '-rf', '%s/%s/star' % (wkdir, self.sample)])
        else:
            os.rename('%s/%s/star/%s.Aligned.sortedByCoord.out.bam' % (wkdir, self.sample, self.sample), '%s/%s/star/%s.bam' % (wkdir, self.sample, self.sample))
    def output(self):
        return luigi.LocalTarget('%s/%s/star/%s.bam' % (wkdir, self.sample, self.sample))


class all_star_align(luigi.Task):
    def requires(self):
        inpt = {}
        for s,p in fastq_dictionary.items():
            inpt[s] = star_align(sample = s, path = p)
        return inpt

    #def run(self):
    #    bam_dict = {x:self.input()[x] for x in fastq_dictionary}
    #    return bam_dict
    def output(self):
        return {x:self.input()[x] for x in fastq_dictionary} 



class featureCounts(luigi.Task):
    sample = luigi.Parameter()
    bam_file = luigi.Parameter()
    

    def requires(self):
        return all_star_align()

    def run(self):
        fc_wkdir = '%s/%s/featureCounts' % (wkdir, self.sample)
        if not os.path.exists(fc_wkdir):
            os.makedirs(fc_wkdir)
        featureCounts_command = [
                                'featureCounts',
                                    '-F', 'GTF',
                                    '-T', '%d' % cores,
                                    '-s', '%d' % stranded,
                                    '-a', '%s' % genome_gtf,
                                    '-o', '%s/%s.prelim.counts' % (fc_wkdir,self.sample),
                                    self.bam_file
                                ]
        fC = subprocess.Popen(featureCounts_command)
        fC.wait()
        if fC.returncode != 0:
            subprocess.call(['rm', '-rf', fc_wkdir])
        else:
            os.rename('%s/%s.prelim.counts' % (fc_wkdir, self.sample), '%s/%s.counts' % (fc_wkdir, self.sample))
            os.rename('%s/%s.prelim.counts.summary' % (fc_wkdir, self.sample), '%s/%s.counts.summary' % (fc_wkdir, self.sample))

    def output(self):
        fc_wkdir = '%s/%s/featureCounts' % (wkdir, self.sample)
        return luigi.LocalTarget('%s/%s.counts' % (fc_wkdir, self.sample))

class all_featureCounts(luigi.Task):

    def requires(self):
        test = all_star_align()
        bam_dict = test.output() 
        return {s:featureCounts(sample = s, bam_file = b.path) for s,b in bam_dict.items()}             

    def run(self):
        #This should give me a dictionary of {sample:gene_counts file}
        return {x:self.input()[x] for x in self.input()}

    def output(self):
        return self.run()

if __name__ == '__main__':
    luigi.run()
