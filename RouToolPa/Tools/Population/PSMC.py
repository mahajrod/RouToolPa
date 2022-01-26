import os
from collections import OrderedDict
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

from RouToolPa.Tools.Abstract import Tool


class PSMC(Tool):
    def __init__(self, path="", max_threads=4, max_memory="100G", max_per_thread_memory="5G"):
        Tool.__init__(self, "plink", path=path, max_threads=max_threads, max_memory=max_memory, max_per_thread_memory=max_per_thread_memory)

    def psmc_plot(self, sample_label_list, psmc_list, generation_time, absolute_mutation_rate,
                  output_prefix, plot_grid=False, min_generations=None, max_generations=None):
        """
        -u FLOAT   absolute mutation rate per nucleotide [2.5e-08]
         -s INT     skip used in data preparation [100]
         -X FLOAT   maximum generations, 0 for auto [0]
         -x FLOAT   minimum generations, 0 for auto [10000]
         -Y FLOAT   maximum popsize, 0 for auto [0]
         -m INT     minimum number of iteration [5]
         -n INT     take n-th iteration (suppress GOF) [20]
         -M titles  multiline mode [null]
         -f STR     font for title, labels and tics [Helvetica,16]
         -g INT     number of years per generation [25]
         -w INT     line width [4]
         -P STR     position of the keys [right top]
         -T STR     figure title [null]
         -N FLOAT   false negative rate [0]
         -S         no scaling
         -L         show the last bin
         -p         convert to PDF (with epstopdf)
         -R         do not remove temporary files
         -G         plot grid
        """

        options = " -g %i" % generation_time
        options += " -u %e" % absolute_mutation_rate
        options += " -p" # TODO: find sence of this options
        options += " -G" if plot_grid else ""
        options += " -M \"%s\"" % ",".join(sample_label_list)
        options += " -X %f" % max_generations if max_generations else ""
        options += " -x %f" % min_generations if min_generations else ""
        options += " %s" % output_prefix
        options += " %s" % " ".join(psmc_list)

        self.execute(options=options, cmd="psmc_plot.pl")

    def prepare_cmd(self, reference_fasta, bam_file, output_prefix, min_coverage, max_coverage, split_dir="split/",
                    min_base_quality=30, min_mapping_quality=30, adjust_mapping_quality=None, min_rms_mapq=30):
        """
        mkdir -p split; time vcfutils.pl splitchr -l 100 ${GENOME}.fai | xargs -I {} -P ${THREADS} sh -c "bcftools mpileup -d 1000000 -q 30 -Q 30 -a AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP -O u -f ${GENOME} -r '{}' ${BAM_LIST}| bcftools call -O u -v -m -f GQ,GP > split/tmp.{}.bcf" && bcftools concat -O u --threads 20 `ls split/tmp.*.bcf | sort -V` | bcftools view -O z -o ${OUTPUT_PREFIX}.vcf.gz - ; rm -r split/
        """
        fq_dir = "%s/fq/" % split_dir

        options = " mkdir -p %s %s;" % (split_dir, fq_dir)
        options += " time vcfutils.pl splitchr -l 100000000000 %s.fai | " % reference_fasta
        options += " xargs -I {} -P %i" % self.threads
        options += " sh -c \"bcftools mpileup "
        options += " -d %i" % (max_coverage * 2)
        options += " -q %i" % min_mapping_quality
        options += " -Q %i" % min_base_quality
        options += " --adjust-MQ %i" % adjust_mapping_quality if adjust_mapping_quality else ""
        options += " -a AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP,SCR,INFO/SCR"
        options += " -O u"
        options += " -f %s " % reference_fasta
        options += " -r '{}' %s| " % bam_file
        options += " bcftools call "
        options += " -c"
        options += " -O u |"
        options += " bcftools view -O v | "
        #options += "" if report_all_positions else " -v"
        #options += " -f GQ |"
        options += "vcfutils.pl vcf2fq -d %i -D %i -Q %i" % (min_coverage, max_coverage, min_rms_mapq)
        options += " > %s/tmp.{}.fq\" && " % fq_dir
        options += "for FILE in `ls %s/tmp.*.fq`; do cat $FILE >> %s.diploid.fq; done" % (fq_dir, output_prefix)

        return options

    def generate_diploid_fastq_from_bam(self, reference_fasta, bam_file, output_prefix,
                                        min_coverage, max_coverage, split_dir="split/",
                                        min_base_quality=30, min_mapping_quality=30,
                                        adjust_mapping_quality=None, min_rms_mapq=30):

        cmd = self.prepare_cmd(reference_fasta, bam_file, output_prefix, min_coverage, max_coverage,
                               split_dir=split_dir, min_base_quality=min_base_quality,
                               min_mapping_quality=min_mapping_quality,
                               adjust_mapping_quality=adjust_mapping_quality, min_rms_mapq=min_rms_mapq)

        self.execute(cmd=cmd, options="")
