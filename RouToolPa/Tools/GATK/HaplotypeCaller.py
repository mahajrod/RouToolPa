#!/usr/bin/env python

__author__ = 'mahajrod'

import os
from RouToolPa.Tools.Abstract import JavaTool
from RouToolPa.Routines import VCFRoutines


class HaplotypeCaller(JavaTool):
    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog=None):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T HaplotypeCaller", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    @staticmethod
    def parse_options_for_parallel_run(reference, alignment, genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
                                       stand_call_conf=30, gvcf_mode=False):

        options = " -R %s" % reference
        options += " -I %s" % alignment
        options += " --genotyping_mode %s" % genotyping_mode if genotyping_mode else ""
        options += " --output_mode %s" % output_mode if output_mode else ""
        #options += " -stand_emit_conf %i" % stand_emit_conf
        options += " -stand_call_conf %i" % stand_call_conf
        options += " --emitRefConfidence GVCF" if gvcf_mode else ""

        return options

    def parse_options(self, reference, alignment, output, genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
                      stand_call_conf=30, gvcf_mode=False, include_region_id_file=None, exclude_region_id_file=None):

        options = " -nct %i" % self.threads
        options += " -R %s" % reference
        options += " -I %s" % alignment
        options += " --genotyping_mode %s" % genotyping_mode if genotyping_mode else ""
        options += " --output_mode %s" % output_mode if output_mode else ""
        #options += " -stand_emit_conf %i" % stand_emit_conf
        options += " -stand_call_conf %i" % stand_call_conf
        options += " --emitRefConfidence GVCF" if gvcf_mode else ""
        options += " -L %s" % include_region_id_file if include_region_id_file else ""
        options += " -XL %s" % exclude_region_id_file if exclude_region_id_file else ""

        options += " -o %s" % output

        return options

    def call(self, reference, alignment, output, genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
             stand_call_conf=30, include_region_id_file=None, exclude_region_id_file=None):
        """
            java -Xmx100g -jar ~/tools/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
              -T HaplotypeCaller \
              -R ${fasta} \
              -I ${bam%bam}realigned.bam \
              --genotyping_mode DISCOVERY \
              --output_mode EMIT_VARIANTS_ONLY \
              -stand_call_conf 30 \
              -o ${bam%bam}raw.vcf
        """
        options = self.parse_options(reference, alignment, output, genotyping_mode=genotyping_mode,
                                     output_mode=output_mode,
                                     stand_call_conf=stand_call_conf, gvcf_mode=False,
                                     include_region_id_file=include_region_id_file,
                                     exclude_region_id_file=exclude_region_id_file)

        self.execute(options)

    def gvcf_call(self, reference, alignment, output, genotyping_mode="DISCOVERY",
                  stand_call_conf=30, include_region_id_file=None, exclude_region_id_file=None):
        """
        java -jar GenomeAnalysisTK.jar \
         -R reference.fasta \
         -T HaplotypeCaller \
         -I sample1.bam \
         --emitRefConfidence GVCF \
         [--dbsnp dbSNP.vcf] \
         [-L targets.interval_list] \
         -o output.raw.snps.indels.g.vcf
        """
        options = self.parse_options(reference, alignment, output, genotyping_mode=genotyping_mode,
                                     stand_call_conf=stand_call_conf, gvcf_mode=True,
                                     include_region_id_file=include_region_id_file,
                                     exclude_region_id_file=exclude_region_id_file)

        self.execute(options)

    def parallel_gvcf_call(self, reference, alignment, output_dir, output_prefix, output,
                           genotyping_mode="DISCOVERY",
                           stand_call_conf=30, max_region_length=1000000, max_seqs_per_region=100,
                           length_dict=None, parsing_mode="parse", region_list=None,
                           region_file_format='simple',
                           remove_intermediate_files=False,
                           gvcf_extension_list=["g.vcf", ],
                           cpus_per_task=1,
                           handling_mode="local",
                           job_name=None,
                           log_prefix=None,
                           error_log_prefix=None,
                           max_running_jobs=None,
                           max_running_time=None,
                           max_memmory_per_cpu=None,
                           modules_list=None,
                           environment_variables_dict=None):
        splited_dir = "%s/splited_gvcf/" % output_dir
        regions_dir = "%s/regions/" % output_dir

        for directory in output_dir, splited_dir:
            self.safe_mkdir(directory)

        region_list, \
            scaffold_to_region_correspondence_dict = self.prepare_region_list_by_length(max_length=max_region_length,
                                                                                        max_seq_number=max_seqs_per_region,
                                                                                        length_dict=length_dict,
                                                                                        reference=None if length_dict is not None else reference,
                                                                                        parsing_mode=parsing_mode,
                                                                                        output_dir=regions_dir,
                                                                                        region_file_format=region_file_format if handling_mode != "slurm" else 'GATK') if region_list is None else region_list

        options = self.parse_options_for_parallel_run(reference, alignment,
                                                      genotyping_mode=genotyping_mode,
                                                      stand_call_conf=stand_call_conf,
                                                      gvcf_mode=True)
        options += " -nct 1"
        options_list = []

        output_index = 1

        output_file_list = []

        if handling_mode == 'local':
            for regions in region_list:
                output_file = "%s/%s_%i.g.vcf" % (splited_dir, output_prefix, output_index)
                region_options = " -o %s" % output_file
                output_file_list.append(output_file)
                #for region in regions:
                #    region_options += " -L %s:%i-%i" % (region[0], region[1], region[2])

                for region in regions:
                    if isinstance(region, str):
                        region_options += " -L %s" % region
                    elif len(region) == 1:
                        region_options += " -L %s" % region[0]
                    elif len(region) == 3:
                        region_options += " -L %s:%i-%i" % (region[0], region[1], region[2])

                options_list.append(options + region_options)
                output_index += 1

            self.parallel_execute(options_list)

            VCFRoutines.combine_same_samples_vcfs(output,
                                                  vcf_list=output_file_list,
                                                  order_vcf_files=True,
                                                  close_fd_after=False,
                                                  extension_list=gvcf_extension_list)

        elif handling_mode == 'slurm':
            number_of_regions = len(region_list)
            region_file = "%s/splited/region_${SLURM_ARRAY_TASK_ID}.list" % regions_dir
            output_file = "%s/%s_${SLURM_ARRAY_TASK_ID}.g.vcf" % (splited_dir, output_prefix)
            options += " -o %s" % output_file
            options += " -L %s" % region_file

            slurm_cmd = "java -Xmx%s -jar %s/%s %s" % (self.max_memory, self.jar_path, self.jar, options)

            last_job_id = self.slurm_run_job(job_name,
                                             log_prefix,
                                             slurm_cmd,
                                             error_log_prefix,
                                             "%s%s.slurm" % (output_dir, output_prefix),
                                             task_index_list=None,
                                             start_task_index=1,
                                             end_task_index=number_of_regions,
                                             max_running_jobs=max_running_jobs,
                                             max_running_time=max_running_time,
                                             cpus_per_task=cpus_per_task,
                                             max_memmory_per_cpu=max_memmory_per_cpu,
                                             modules_list=modules_list,
                                             environment_variables_dict=environment_variables_dict)

            print("Submitted job  %s" % last_job_id)

            #VCFRoutines.combine_same_samples_vcfs(output_file_list,
            #                                      output,
            #                                      close_fd_after=False,
            #                                      extension_list=gvcf_extension_list)
        else:
            print("ERROR!!! Unrecognized handling mode!")
