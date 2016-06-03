import luigi
import basic_rna_pipe as brp
import eisa_rna_pipe as erp
import other_features_rna_pipe as ofrp

class all_some_task(luigi.WrapperTask):
    require = luigi.TaskParameter()

    def requires(self):
        yield {sample: self.require(sample=sample)
               for sample in brp.fastqs().output()}
