import luigi
import subprocess
import json
import os

'''
The two functions below are used to convert between paths on the host system
and the path in the docker data volume.

These functions will not be needed once luigi itself is put into docker
'''


def host_path(data_volume, docker_path):
    '''Used to map relative paths to absolute paths in docker data volume
       Allows luigi's localTarget to register
    '''
    cmd = ['docker', 'inspect', data_volume]
    # Convert json output to string
    inspect_string = subprocess.check_output(cmd).decode('utf-8')
    # Convert to json
    inspect_json = json.loads(inspect_string)
    # Retrieve host path from data volume info
    host_path_suffix = inspect_json[0]['Mounts'][0]['Source']
    return os.path.join(host_path_suffix, docker_path)


def docker_path(data_volume, real_path):
    cmd = ['docker', 'inspect', data_volume]
    # Convert json output to string
    inspect_string = subprocess.check_output(cmd).decode('utf-8')
    # Convert to json
    inspect_json = json.loads(inspect_string)
    # Retrieve mount point in data volume
    dest_mount = inspect_json[0]['Mounts'][0]['Destination']
    # Get docker part of real path
    docker_p = real_path.rsplit("_data/")[-1]
    return os.path.join(dest_mount, docker_p)


class globalConfig(luigi.Config):
    experiment_name = luigi.Parameter(default='test')
    data_volume = luigi.Parameter(default='rna_datastore')
    paired_end = luigi.BoolParamter()
    fastq_dict = luigi.DictParameter()


class kallisto_index(luigi.Task):
    transcript_fasta = luigi.Parameter()
    index_name = transcript_fasta.split(".fa")[0] + '.idx'

    def run(self):
        docker_cmd = ['docker', 'run',
                      '--volumes-from', globalConfig().data_volume,
                      'kallisto', 'index',
                      '-i', '/rna_data/index/%s' % self.index_name,
                      '/rna_data/transcript_fastas/%s'
                      % globalConfig().transcript_fasta]

        subprocess.call(docker_cmd)

    def output(self):
        h_path = host_path(data_volume=globalConfig().data_volume,
                           volume_path='index/%s'
                           % self.index_name)
        return luigi.LocalTarget(h_path)


class kallisto_quant(luigi.Task):
    sample_fastq = luigi.Parameter()
    paired_end = luigi.Parameter()

    def requires(self):
        return kallisto_index

    def run(self):
        if self.paried_end:
            docker_cmd = ['docker', 'run',
                          '--volumes-from', globalConfig().data_volume,
                          'kallisto', 'quant',
                          '-i', docker_path(self.input().path),
                          '-o', '/rna_data/%s'
                          % globalConfig().experiment_name,
                          globalConfig().fastq_dict[self.sample]
                          ]
        else:
            docker_cmd = ['docker', 'run',
                          '--volumes-from', globalConfig().data_volume,
                          'kallisto', 'quant',
                          '-i', docker_path(self.input().path),
                          '-o', '/rna_data/%s'
                          % globalConfig().experiment_name,
                          '--single', '-l', '200', '-s', '20',
                          self.sample_fastq]

        subprocess.call(docker_cmd)

    def output():
        h_path = host_path(data_volume=globalConfig().data_volume,
                           volume_path=globalConfig().experiment_name)
        return luigi.LocalTarget(h_path)
