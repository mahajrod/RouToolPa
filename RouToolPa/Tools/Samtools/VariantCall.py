#!/usr/bin/env python


from RouToolPa.Tools.Abstract import Tool


class VariantCall(Tool):
    """
    Class for samtools 1.0+
    Several subcommands are not implemented
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "samtools", path=path, max_threads=max_threads)

        # bam/sam flag values:
        self.bam_flags = {
                          "read_paired": 1,
                          "read_mapped_in_proper_pair": 2,
                          "read_unmapped": 4,
                          "mate_unmapped": 8,
                          "read_reverse_strand": 16,
                          "mate_reverse_strand": 32,
                          "first_in_pair": 64,
                          "second_in_pair": 128,
                          "not_primary_alignment": 256,
                          "read_fails_platform/vendor_quality_checks": 512,
                          "read_is_PCR_or_optical_duplicate": 1024,
                          "supplementary_alignment": 2048
                          }

    def prepare_cmd(self, reference_fasta, bam_list, output_prefix, chunk_length=1000000, split_dir="split/", max_coverage=None,
                    min_base_quality=30, min_mapping_quality=30):
        """
        mkdir -p split; time vcfutils.pl splitchr -l 100 ${GENOME}.fai | xargs -I {} -P ${THREADS} sh -c "bcftools mpileup -d 1000000 -q 30 -Q 30 -a AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP -O u -f ${GENOME} -r '{}' ${BAM_LIST}| bcftools call -O u -v -m -f GQ,GP > split/tmp.{}.bcf" && bcftools concat -O u --threads 20 `ls split/tmp.*.bcf | sort -V` | bcftools view -O z -o ${OUTPUT_PREFIX}.vcf.gz - ; rm -r split/
        """

        options = " mkdir -p %s;" % split_dir
        options += " time vcfutils.pl splitchr -l %i %s.fai | " % (chunk_length, reference_fasta)
        options += " xargs -I {} -P %i" % self.threads
        options += " sh -c \"bcftools mpileup "
        options += " -d %i" % max_coverage if max_coverage else ""
        options += " -q %i" % min_mapping_quality
        options += " -Q %i" % min_base_quality
        options += " -a AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP"
        options += " -O u"
        options += " -f %s " % reference_fasta
        options += " -r '{}' %s| " % " ".join(bam_list)
        options += " bcftools call -O u -v -m -f GQ,GP > split/tmp.{}.bcf\" &&"
        options += " bcftools concat -O u --threads %i `ls %s/tmp.*.bcf | sort -V` |" % (self.threads, split_dir)
        options += " bcftools view -O z -o %s.vcf.gz -; " % output_prefix
        options += " rm -r %s" % split_dir

        return options

    def call_variants(self, reference_fasta, output_prefix, bam_list, chunk_length=1000000, split_dir="split/", max_coverage=None,
                      min_base_quality=30, min_mapping_quality=30):

        cmd = self.prepare_cmd(reference_fasta, bam_list, output_prefix, chunk_length=chunk_length,
                               split_dir=split_dir, max_coverage=max_coverage,
                               min_base_quality=min_base_quality, min_mapping_quality=min_mapping_quality)

        self.execute(options="", cmd=cmd)
