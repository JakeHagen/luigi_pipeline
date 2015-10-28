import luigi
import subprocess
import os
import rpy2.robjects as robjects                                            
import pandas as pd 

class parameters(luigi.Task):
    """class to take all parameters
    fastq_file needs to be in the format specified below, if the experiment is 
    paired end, the second fastq file should be after the first, only seperated 
    by a space
    
    s01:/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi1_CAGATC_L008_R1_001.C6673ACXX.fastq.gz
    s02:/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi2_ATCACG_L008_R1_001.C6673ACXX.fastq.gz
    s03:/sc/orga/projects/argmac01a/QC_C211.B857_huh7_DCPS.SE.RNASeqRibozero.RAPiD.Human/fastqs/lsi3_TCGGCA_L008_R1_001.C6673ACXX.fastq.gz
    """
    
    fastq_file = luigi.Parameter(
                    default = '/hpc/users/hagenj02/luigi_pipeline/fastqfile')
    
    wkdir = luigi.Parameter(default = os.getcwd())

    genome_fasta = luigi.Parameter(
            default = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/'
                        'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
                                   )
    genome_gtf = luigi.Parameter(
            default = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/'
                        'Homo_sapiens.GRCh37.75.gtf'
                                 )
    star_genome_folder = luigi.Parameter(
            default = '/sc/orga/projects/Houton_Sander/genomes/ensembl37/'
                        'star_genome_ensembl37_99'
                                         )
    read_length = luigi.IntParameter(default = 100)
    stranded = luigi.IntParameter(default = 0)
    paried = luigi.IntParameter(default = 0)                                            
    cores = luigi.IntParameter(default = 6)    


class fastqs(luigi.Task):
    
    def requires(self):
        return parameters()
    
    def run(self):    
        fastq_dict = {}
        with open(parameters().fastq_file) as f:
            for line in f:
                sample,path = line.strip().split(":")
                fastq_dict[sample] = path
        return fastq_dict
    
    def output(self):
        fastq_dict = self.run()
        return {sample:luigi.LocalTarget(path) for sample,path in fastq_dict.items()} 


class index_STAR_genome(luigi.Task):
    
    def requires(self):
        return parameters()

    def run(self):
        if not os.path.exists(parameters().star_genome_folder):
            os.makedirs(parameters().star_genome_folder) 
        star_command = [
                        'STAR',
                        '--runThreadN %d' % parameters().cores,
                        '--runMode genomeGenerate',
                        '--genomeDir %s' % parameters().star_genome_folder,
                        '--genomeFastaFiles %s' % parameters().genome_fasta,
                        '--sjdbGTFfile %s' % parameters().genome_gtf,
                        '--sjdbOverhang %d' % (parameters.read_length - 1)
                        ]
        star = subprocess.Popen(star_command)
        star.wait()
        if star.returncode != 0:
            subprocess.call(['rm -rf %s' % parameters().star_genome_folder])
        os.rename('Log.out', '%s/Log.out' % parameters().star_genome_folder)
    def output(self):
        return luigi.LocalTarget('%s/Log.out' % parameters().star_genome_folder)


class wk_dir(luigi.Task):
    
    def requires(self):
        return index_STAR_genome()
    
    def run(self):
        os.makedirs(parameters().wkdir)
    
    def output(self):
        return luigi.LocalTarget(parameters().wkdir)


class star_align(luigi.Task):
    
    sample = luigi.Parameter()
    file_location = luigi.Parameter()
    def requires(self):
        return index_STAR_genome(), wk_dir(), fastqs()

    def run(self):
        if not os.path.exists('%s/%s/star' % (parameters().wkdir, self.sample)):                    os.makedirs('%s/%s/star' % (parameters().wkdir, self.sample))
        s_command = [
                'STAR',
                    '--genomeDir %s' % parameters().star_genome_folder,
                    #'--sjdbGTFfile %s' % genome_gtf,
                    '--readFilesIn %s' % self.file_location.path,
                    '--readFilesCommand zcat',
                    '--runThreadN %d' % parameters().cores,
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
                    '--outFileNamePrefix %s/%s/star/%s.' % 
                        (parameters().wkdir, self.sample, self.sample)
                    ]
        star = subprocess.Popen(s_command)
        star.wait()

        '''Below removes folder if star command was unsuccessful 
        (files remain when star fails)
        renaming bam insures that if batch job gets cut while star is running
        the task will not be marked as complete 
        '''
        if star.returncode != 0:
            subprocess.call(['rm', '-rf', '%s/%s/star' 
                                % (parameters().wkdir, self.sample)])
        else:
            os.rename('%s/%s/star/%s.Aligned.sortedByCoord.out.bam' 
                        % (parameters().wkdir, self.sample, self.sample), 
                        '%s/%s/star/%s.bam' 
                            % (parameters().wkdir, self.sample, self.sample))
    def output(self):
        return luigi.LocalTarget('%s/%s/star/%s.bam' 
                                    % (parameters().wkdir, self.sample, self.sample))


class all_star_align(luigi.Task):
    '''Requires splits out samples and runs star align concurrently
    output is dictionary with sample name and luigi target of bam file
    '''
    def requires(self):
        fastq_dict = fastqs().output()
        return {
                s:star_align(sample = s, file_location = p) 
                    for s,p in fastqs_dict.items()
                }
    def output(self):       
        bam_dict = self.input()
        return bam_dict 


class featureCounts(luigi.Task):
    sample = luigi.Parameter()
    bam_file = luigi.Parameter()
    
    def requires(self):
        return all_star_align()

    def run(self):
        fc_wkdir = '%s/%s/featureCounts' % (parameters().wkdir, self.sample)
        if not os.path.exists(fc_wkdir):
            os.makedirs(fc_wkdir)
        featureCounts_command = [
                                'featureCounts',
                                    '-F', 'GTF',
                                    '-T', '%d' % cores,
                                    '-s', '%d' % stranded,
                                    '-a', '%s' % genome_gtf,
                                    '-o', '%s/%s.prelim.counts' % (fc_wkdir,self.sample),
                                    self.bam_file.path
                                ]
        fC = subprocess.Popen(featureCounts_command)
        fC.wait()
        if fC.returncode != 0:
            subprocess.call(['rm', '-rf', fc_wkdir])
        else:
            os.rename('%s/%s.prelim.counts' % (fc_wkdir, self.sample), 
                        '%s/%s.counts' % (fc_wkdir, self.sample))
            os.rename('%s/%s.prelim.counts.summary' % (fc_wkdir, self.sample), 
                        '%s/%s.counts.summary' % (fc_wkdir, self.sample))

    def output(self):
        fc_wkdir = '%s/%s/featureCounts' % (parameters().wkdir, self.sample)
        return luigi.LocalTarget('%s/%s.counts' % (fc_wkdir, self.sample))


class all_featureCounts(luigi.Task):
    '''requires splits out samples and run featureCounts concurrently
    output returns dictionary with sample name and count file
    '''
    def requires(self):
        bam_dict = all_star_align().output() 
        return {
                s:featureCounts(sample = s, bam_file = b.path) 
                    for s,b in bam_dict.items()
                }             

    def output(self):
        counts_dict = self.input()
        return counts_dict


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
       # pandas2ri.activate()
       # r = robjects.r
       # robjects.globalenv['experimentGroups'] = robjects.StrVector(experiment_group)
       
       #robjects.globalenv['countTable'] = pandas2ri.py2ri(count_table)
        
        limma = robjects.r('''
                library("limma")
                library("edgeR")
                experimentGroups <- c("c13", "c14", "c15", "s01", "s02", "s03")
                design <- model.matrix(~0+factor(experimentGroups))
                colnames(design) <- c(unique(experimentGroups))
                
                countTable <- read.csv("counts")
                rownames(countTable) <- countTable$Gene
                countTable <- countTable[c(2,3,4,5,6,7)]
                t <- DGEList(countTable, group = experimentGroups)

                keep <- rowSums(cpm(t)>1) >= 3  
                y <- t[keep,] #throw out genes that do not have at least one cpm in at least 3 samples 

                dge <- calcNormFactors(y) #Use TMM normalization to correct for library size
                v <- voom(dge, design=design,plot=FALSE) #Voom transformation, convert to log2CPM and associate a weight based on variance

                fit <- lmFit(v,design = design) #Fit glm model
                contrast.matrix <- makeContrasts(siRNA-control, levels=design)#Contrast matrix for constrasts of interest
                fit2 <- contrasts.fit(fit,contrast.matrix)
                fit2 <- eBayes(fit2)

                # Output list of all genes for each contrast that had counts
                # Adjusted p values by BH, (BH is default)
                siRNA_control_full <- topTable(fit2, number = 200000, coef=1)
                siRNA_control_full$gene <- rownames(siRNA_control_full)
                siRNA_control_padj.05 <- topTable(fit2, number = 200000, coef=1, adjust="BH", p.value = .1)
                write.table(siRNA_control_padj.05, file = "/hpc/users/hagenj02/luigi_pipeline/deg_test.txt", col.names = NA, sep = "\t", quote = FALSE)
            ''')
        limma()
        #r['source']("/hpc/users/hagenj02/luigi_pipeline/script.R")
    def output(self):
        return luigi.LocalTarget('/hpc/users/hagenj02/luigi_pipeline/deg_test.txt')
        
    #def output():
    #    return [self.input()[y].path for y in self.input()]

if __name__ == '__main__':
    luigi.run()
