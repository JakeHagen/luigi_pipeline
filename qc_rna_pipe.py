import luigi
import os
import subprocess
import basic_rna_pipe as brp


class fastqc(luigi.Task):
    '''Runs fastqc, quality control program on fastqs
    '''
    sample = luigi.Parameter()

    def requires(self):
        return brp.fastqs(sample=self.sample)

    def output_dir(self):
        return '%s/fastqc' % brp.configs().exp_dir

    def run(self):
        try:
            os.makedirs(self.output_dir())
        except OSError:
            pass

        fastqc_command = ['fastqc', self.input().path,
                          '-o', self.output_dir()]
        subprocess.call(fastqc_command)

    def output(self):
        file_name = os.path.basename(self.input().path).split(".")[0]
        return luigi.LocalTarget('%s/%s_fastqc.html' %
                                (self.output_dir(), file_name))
                                
