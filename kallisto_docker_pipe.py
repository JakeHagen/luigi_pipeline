import luigi
import subprocess
import json
import os

experiment_name = 'experiment_test'
transcript_fasta = 'gencode.vM9.transcripts.fa.gz'
data_volume = 'rna_datastore'


def docker_real_path(data_volume, volume_path):
    '''Used to map relative paths to absolute paths in docker data volume
       Allows luigi's localTarget to register
    '''
    cmd = ['docker', 'inspect', data_volume]
    # Convert json output to string
    inspect_string = subprocess.check_output(cmd).decode('utf-8')
    # Convert to json
    inspect_json = json.loads(inspect_string)
    # Retrieve host path from data volume info
    host_path = inspect_json[0]['Mounts'][0]['Source']
    return os.path.join(host_path, volume_path)


class kallisto_index(luigi.Task):
    index_name = transcript_fasta.split(".fa")[0] + '.idx'

    def run(self):
        docker_cmd = ['docker', 'run',
                      '--volumes-from', data_volume,
                      'kallisto', 'index',
                      '-i', '/rna_data/index/%s' % self.index_name,
                      '/rna_data/transcript_fastas/%s' % transcript_fasta]

        subprocess.call(docker_cmd)

    def output(self):
        abs_path = docker_real_path(data_volume=data_volume,
                                    volume_path='rna_data/index/%s'
                                    % self.index_name)
        return luigi.LocalTarget(abs_path)


class kallisto_quant(luigi.Task):
    sample_fastq = luigi.Parameter()

    def requires(self):
        return kallisto_index

    def run(self):
        docker_cmd = ['docker', 'run',
                      
                      'kallisto', 'quant',
                      '-i', self.input().path,
                      '-o', '%s/luigi_out' % work_dir,
                      '--single', '-l', '200', '-s', '20',
