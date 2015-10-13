import luigi
import subprocess
import os
import readline
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import pandas as pd

#class params(luigi.Task):

    #will output will be, can be pointed at an existing folder or the pipeline will create it
wkdir = '/sc/orga/projects/argmac01a/luigi_rna_ensembl37'
    # Only needed if STAR genome needs to be generated
genome_fasta = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa' 
genome_gtf = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/Homo_sapiens.GRCh37.75.gtf'

"""Either where the STAR genome is located or were it will be created
IMPORTANT: if the STAR genome needs to be created, do not create the folder,
Let the pipeline create it."""
star_genome_folder = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/star_genome_ensembl37_99'

    # RNAseq experiment info
read_length = 100
stranded = 0
    #paried = ?
    #prep = polyA or ribozero

    # Should be multiple of samples
cores = 6

"""
s01:/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi1_CAGATC_L008_R1_001.C6673ACXX.fastq.gz
s02:/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi2_ATCACG_L008_R1_001.C6673ACXX.fastq.gz
s03:/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi3_TCGGCA_L008_R1_001.C6673ACXX.fastq.gz
"""

class params(luigi.Task):
    
    fq = luigi.Parameter()
    wkdir = luigi.Parameter(default = '/sc/orga/projects/argmac01a/luigi_rna_ensembl37')                       
    # Only needed if STAR genome needs to be generated                          
    genome_fasta = luigi.Parameter(
            default = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
                                   )
    genome_gtf = luigi.Parameter(
            default = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/Homo_sapiens.GRCh37.75.gtf'
                                 )
    star_genome_folder = luigi.Parameter(
            default = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/star_genome_ensembl37_99'
                                         )
    read_length = luigi.IntParameter(default = 100)
    stranded = luigi.IntParameter(default = 0)
    #paried = ?
    #prep = polyA or ribozero                                             
    cores = luigi.IntParameter(default = 6)    

pm = params

class fastqs(luigi.Task):
    
    def requires(self):
        return params()
    def output(self):
        d = {}
        with open(pm().fq) as f:
            for line in f:
                key,val = line.strip().split(":")
                d[key] = val  
        return d

class index_STAR_genome(luigi.Task):
    
    def requires(self):
        return fastqs()

    def run(self):
        if not os.path.exists(pm().star_genome_folder):
            os.makedirs(pm().star_genome_folder) 
        star_command = [
                        'STAR',
                        '--runThreadN %d' % pm.cores,
                        '--runMode genomeGenerate',
                        '--genomeDir %s' % pm().star_genome_folder,
                        '--genomeFastaFiles %s' % pm().genome_fasta,
                        '--sjdbGTFfile %s' % pm().genome_gtf,
                        '--sjdbOverhang %d' % (pm.read_length - 1)
                        ]
        star = subprocess.Popen(star_command)
        star.wait()
        if star.returncode != 0:
            subprocess.call(['rm -rf %s' % pm().star_genome_folder])
        os.rename('Log.out', '%s/Log.out' % pm().star_genome_folder)
    def output(self):
        return luigi.LocalTarget('%s/Log.out' % pm().star_genome_folder)


class wk_dir(luigi.Task):
    
    def requires(self):
        return index_STAR_genome()
    def run(self):
        os.makedirs(pm().wkdir)
    def output(self):
        return luigi.LocalTarget(pm().wkdir)


class star_align(luigi.Task):
    
    sample = luigi.Parameter()
    path = luigi.Parameter()
    def requires(self):
        return index_STAR_genome(), wk_dir()

    def run(self):
        os.makedirs('%s/%s/star' % (pm().wkdir, self.sample))
        s_command = [
                'STAR',
                    '--genomeDir %s' % pm().star_genome_folder,
                    #'--sjdbGTFfile %s' % genome_gtf,
                    '--readFilesIn %s' % self.path,
                    '--readFilesCommand zcat',
                    '--runThreadN %d' % pm.cores, #(cores/len(fastq_dictionary)),
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
                    '--outFileNamePrefix %s/%s/star/%s.' % (pm().wkdir, self.sample, self.sample)
                    ]
        star = subprocess.Popen(s_command)
        star.wait()
        # Below removes folder if star command was unsuccessful 
        # (files remain when star fails)
        # renaming bam insures that if batch job gets cut while star is running
        # the task will not be marked as complete 
        if star.returncode != 0:
            subprocess.call(['rm', '-rf', '%s/%s/star' % (pm().wkdir, self.sample)])
        else:
            os.rename('%s/%s/star/%s.Aligned.sortedByCoord.out.bam' % (pm().wkdir, self.sample, self.sample), '%s/%s/star/%s.bam' % (pm().wkdir, self.sample, self.sample))
    def output(self):
        return luigi.LocalTarget('%s/%s/star/%s.bam' % (pm().wkdir, self.sample, self.sample))


class all_star_align(luigi.Task):
    
    def requires(self):
        inpt = {}
        for s,p in fastqs().output().items():
            inpt[s] = star_align(sample = s, path = p)
        return inpt

    def output(self):       
        return {x:self.input()[x] for x in self.input()} 


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
        bam_dict = all_star_align().output() 
        return {s:featureCounts(sample = s, bam_file = b.path) for s,b in bam_dict.items()}             

    def output(self):        
        return {x:self.input()[x] for x in self.input()}


class diff_exp_analysis(luigi.Task):
    
    def requires(self):
        return all_featureCounts()

    def run(self):
        sample_names = [x for x in self.input()]
        count_files = [self.input()[y].path for y in self.input()]
        experiment_group = [x[0] for x in self.input()]
        
        files = [
                pd.read_table(self.input()[name].path, 
                                skiprows=2, 
                                index_col=0,
                                names = ['Gene', 't', 'e', 's', 'r', 'w', name],
                                usecols = ['Gene', name], 
                                header=None)
                    for name in self.input()
                    ]
        count_table = pd.concat(files, axis = 1).sort_index(axis=1)
       # count_table.to_csv("/hpc/users/hagenj02/luigi_pipeline/counts")
        pandas2ri.activate()
        r = robjects.r
        robjects.globalenv["experimentGroups"] = robjects.StrVector(experiment_group)
        robjects.globalenv["countTable"] = pandas2ri.py2ri(count_table)
        
        r['source']("/hpc/users/hagenj02/luigi_pipeline/script.R")
    def output(self):
        return luigi.LocalTarget('/hpc/users/hagenj02/luigi_pipeline/deg_test.txt')
        
    #def output():
    #    return [self.input()[y].path for y in self.input()]

if __name__ == '__main__':
    luigi.run()
