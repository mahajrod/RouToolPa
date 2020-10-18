#!/usr/bin/env python

__author__ = 'mahajrod'
import os
import shutil
from RouToolPa.Collections.General import IdList
from RouToolPa.Tools.Abstract import Tool

from RouToolPa.Routines import VCFRoutines


class HaplotypeCaller4(Tool):
    def __init__(self,  max_threads=4,  max_memory=None, timelog=None):
        Tool.__init__(self,
                      "gatk HaplotypeCaller",
                      max_threads=max_threads, max_memory=max_memory, timelog=timelog)

    @staticmethod
    def parse_options_for_parallel_run(reference, alignment, output_mode=None,
                                       stand_call_conf=30, gvcf_mode=False,
                                       ignore_softclipped_bases=False):

        options = " -R %s" % reference

        if isinstance(alignment, str):
            options += " -I %s" % alignment
        else:
            for alignment_file in alignment:
                options += " -I %s" % alignment_file

        options += " --output_mode %s" % output_mode if output_mode else ""
        #options += " -stand_emit_conf %i" % stand_emit_conf
        options += " --standard-min-confidence-threshold-for-calling %i" % stand_call_conf
        options += " --emit-ref-confidence GVCF" if gvcf_mode else ""
        options += " --dont-use-soft-clipped-bases" if ignore_softclipped_bases else ""

        return options

    def parse_options(self, reference, alignment, output, output_mode=None,
                      stand_call_conf=30, gvcf_mode=False, include_region_id_file=None, exclude_region_id_file=None):

        #options = " -nct %i" % self.threads
        options = " -R %s" % reference
        if isinstance(alignment, str):
            options += " -I %s" % alignment
        else:
            for alignment_file in alignment:
                options += " -I %s" % alignment_file
        options += " --output_mode %s" % output_mode if output_mode else ""
        #options += " -stand_emit_conf %i" % stand_emit_conf
        options += " --standard-min-confidence-threshold-for-calling %i" % stand_call_conf
        options += " --emit-ref-confidence GVCF" if gvcf_mode else ""
        options += " -L %s" % include_region_id_file if include_region_id_file else ""
        options += " -XL %s" % exclude_region_id_file if exclude_region_id_file else ""

        options += " -O %s" % output

        return options

    def call(self, reference, alignment, output, output_mode=None,
             stand_call_conf=30, include_region_id_file=None, exclude_region_id_file=None):
        """
            gatk HaplotypeCaller \
              -R ${fasta} \
              -I ${bam%bam}realigned.bam \
              --output_mode EMIT_VARIANTS_ONLY \
              -stand_call_conf 30 \
              -o ${bam%bam}raw.vcf
        """
        options = self.parse_options(reference, alignment, output,
                                     output_mode=output_mode,
                                     stand_call_conf=stand_call_conf, gvcf_mode=False,
                                     include_region_id_file=include_region_id_file,
                                     exclude_region_id_file=exclude_region_id_file)

        self.execute(options,
                     cmd=("gatk --java-options -Xmx%s HaplotypeCaller" % self.max_memory) if self.max_memory else None)

    def gvcf_call(self, reference, alignment, output,
                  stand_call_conf=30, include_region_id_file=None, exclude_region_id_file=None):
        """
        gatk HaplotypeCaller \
         -R reference.fasta
         -I sample1.bam \
         --emit-ref-confidence GVCF \
         [--dbsnp dbSNP.vcf] \
         [-L targets.interval_list] \
         -o output.raw.snps.indels.g.vcf
        """
        options = self.parse_options(reference, alignment, output,
                                     stand_call_conf=stand_call_conf, gvcf_mode=True,
                                     include_region_id_file=include_region_id_file,
                                     exclude_region_id_file=exclude_region_id_file)

        self.execute(options,
                     cmd=("gatk --java-options -Xmx%s HaplotypeCaller" % self.max_memory) if self.max_memory else None)

    def parallel_call(self, reference, alignment, output_dir, output_prefix,
                      stand_call_conf=30, max_region_length=1000000, max_seqs_per_region=100,
                      length_dict=None, parsing_mode="parse", region_list=None,
                      region_file_format='simple',
                      remove_intermediate_files=False,
                      cpus_per_task=1,
                      handling_mode="local",
                      job_name=None,
                      log_prefix=None,
                      error_log_prefix=None,
                      max_running_jobs=None,
                      max_running_time=None,
                      max_memmory_per_cpu=None,
                      modules_list=None,
                      environment_variables_dict=None,
                      black_list_scaffold_id_file=None,
                      gvcf_mode=False,
                      ignore_softclipped_bases=False):

        splited_dir = "%s/splited/" % output_dir
        regions_dir = "%s/regions/" % output_dir

        from RouToolPa.Tools.GATK4 import SortVcf4
        sequence_dict = reference[:-5] + "dict"
        SortVcf4.max_memory = self.max_memory
        SortVcf4.path = self.path

        for directory in output_dir, splited_dir:
            self.safe_mkdir(directory)

        if black_list_scaffold_id_file:
            if isinstance(black_list_scaffold_id_file, str):
                black_scaffolds_list = IdList(filename=black_list_scaffold_id_file)
            else:
                black_scaffolds_list = black_list_scaffold_id_file
        else:
            black_scaffolds_list = []
        region_list, \
            scaffold_to_region_correspondence_dict = self.prepare_region_list_by_length(max_length=max_region_length,
                                                                                        max_seq_number=max_seqs_per_region,
                                                                                        length_dict=length_dict,
                                                                                        reference=None if length_dict is not None else reference,
                                                                                        parsing_mode=parsing_mode,
                                                                                        output_dir=regions_dir,
                                                                                        black_list_scaffolds=black_scaffolds_list,
                                                                                        region_file_format=region_file_format if handling_mode != "slurm" else 'GATK') if region_list is None else region_list

        options = self.parse_options_for_parallel_run(reference, alignment,

                                                      stand_call_conf=stand_call_conf,
                                                      gvcf_mode=gvcf_mode,
                                                      ignore_softclipped_bases=False)
        #options += " -nct 1"
        options_list = []

        output_index = 1

        output_file_list = []

        output_extension = "g.vcf" if gvcf_mode else "vcf"

        if handling_mode == 'local':
            for regions in region_list:

                output_file = "%s/%s_%i.%s" % (splited_dir, output_prefix, output_index, output_extension)
                region_options = " -O %s" % output_file
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
            print("Variant calling....")
            self.parallel_execute(options_list,
                                  cmd=("gatk --java-options -Xmx%s HaplotypeCaller" % self.max_memory) if self.max_memory else None)
            unsorted_combined_vcf = "%s/%s.unsorted.%s" % (output_dir, output_prefix, output_extension)
            sorted_combined_vcf = "%s/%s.%s" % (output_dir, output_prefix, output_extension)
            print("Combining variants...")
            VCFRoutines.combine_same_samples_vcfs(unsorted_combined_vcf,
                                                  vcf_list=output_file_list,
                                                  order_vcf_files=True,
                                                  close_fd_after=False,
                                                  extension_list=[".vcf",])
            print("Sorting...")
            SortVcf4.sort_vcf(unsorted_combined_vcf, sorted_combined_vcf, sequence_dict)
            shutil.rmtree(splited_dir)
            shutil.rmtree(regions_dir)
            os.remove(unsorted_combined_vcf)

        elif handling_mode == 'slurm':
            number_of_regions = len(region_list)
            region_file = "%s/splited/region_${SLURM_ARRAY_TASK_ID}.list" % regions_dir
            output_file = "%s/%s_${SLURM_ARRAY_TASK_ID}.%s" % (splited_dir, output_prefix, output_extension)
            options += " -O %s" % output_file
            options += " -L %s" % region_file

            slurm_cmd = "gatk --java-options -Xmx%s HaplotypeCaller" % self.max_memory if self.max_memory else "gatk HaplotypeCaller"
            slurm_cmd += " %s" % options

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
