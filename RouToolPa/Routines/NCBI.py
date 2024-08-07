#!/usr/bin/env python
import os
import re
import sys
import time
from collections.abc import Iterable
from collections import OrderedDict

if sys.version_info[0] == 2:
    from urllib2 import URLError
else:
    from urllib.error import URLError

import xmltodict
import numpy as np
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from RouToolPa.Routines.Sequence import SequenceRoutines
from RouToolPa.Routines.Annotations import AnnotationsRoutines
from RouToolPa.Parsers.GFF import CollectionGFF
from RouToolPa.Collections.General import IdList, SynDict, TwoLvlDict, IdSet


class AssemblySummary(OrderedDict):

    def __init__(self, assembly_summary_biopython_dict_element):
        OrderedDict.__init__(self)
        #parameters_dict = entrez_summary_biopython['DocumentSummarySet']['DocumentSummary'][0]
        self.stats = self.parse_stats_from_assembly_summary_metadata(assembly_summary_biopython_dict_element['Meta'])

        for parameter in 'SpeciesName', 'Organism', 'Taxid', 'AssemblyName', 'FtpPath_GenBank', 'PropertyList', \
                         'SubmissionDate', 'LastUpdateDate', 'AssemblyAccession', 'LastMajorReleaseAccession', \
                         'AssemblyType', 'PartialGenomeRepresentation', 'AssemblyClass', 'AnomalousList',\
                         'AssemblyStatus', 'BioSampleId', 'BioSampleAccn', 'WGS', 'Biosource':
            if parameter in assembly_summary_biopython_dict_element:
                setattr(self, parameter, assembly_summary_biopython_dict_element[parameter])
            else:
                setattr(self, parameter, None)
                #exec("self[%s] = None" % parameter)

        if self.Biosource:
            #print self.Biosource
            self['Isolate'] = self.Biosource['Isolate'] if self.Biosource['Isolate'] else None
            if ('InfraspeciesList' in self.Biosource) and self.Biosource['InfraspeciesList']:
                self['Sub_type'] = self.Biosource['InfraspeciesList'][0]['Sub_type'] if 'Sub_type' in self.Biosource['InfraspeciesList'][0] else None
                self['Sub_value'] = self.Biosource['InfraspeciesList'][0]['Sub_value'] if 'Sub_value' in self.Biosource['InfraspeciesList'][0] else None
            else:
                self['Sub_type'] = None
                self['Sub_value'] = None
            if ('Sex' in self.Biosource) and self.Biosource['Sex']:
                self['Sex'] = self.Biosource['Sex']
            else:
                self['Sex'] = None

        else:
            self['Isolate'] = None
            self['Sub_type'] = None
            self['Sub_value'] = None
            self['Sex'] = None

        for parameter in 'Isolate', 'Sub_type', 'Sub_value':
            setattr(self, parameter, self[parameter])

        for entry in self.stats:
            self[entry] = self.stats[entry]
        for entry in assembly_summary_biopython_dict_element:
            self[entry] = assembly_summary_biopython_dict_element[entry]

        for parameter in "alt_loci_count", "chromosome_count", "contig_count", "contig_l50", "contig_n50", \
                         "non_chromosome_replicon_count", "replicon_count", "scaffold_count_all", \
                         "scaffold_count_placed", "scaffold_count_unlocalized", "scaffold_count_unplaced",\
                         "scaffold_l50", "scaffold_n50", "total_length", "ungapped_length":
            setattr(self, parameter, self[parameter])

        self.full_genome = True if 'full-genome-representation' in self.PropertyList else False
        self.chromosome_lvl = True if 'has-chromosome' in self.PropertyList else False
        self.chloroplast = True if 'has-chloroplast' in self.PropertyList else False
        self.mitochondrion = True if 'has-mitochondrion' in self.PropertyList else False
        self.plasmid = True if 'has-plasmid' in self.PropertyList else False

        self.header_list = ['SpeciesName', 'Taxid', 'Isolate', 'Sub_type', 'Sub_value', 'Sex', 'BioSampleId', 'BioSampleAccn',
                            'WGS', 'AssemblyAccession', 'AssemblyName', 'SubmissionDate', 'LastUpdateDate',
                            'AssemblyType', 'AssemblyStatus',
                            'WholeGenome', "ChromosomeLvl", 'Chloroplast', "Mitochondrion", "Plasmid",
                            'total_length', 'ungapped_length', "alt_loci_count",
                            'chromosome_count', 'non_chromosome_replicon_count', 'replicon_count',
                            'scaffold_count_all', 'scaffold_count_placed', 'scaffold_count_unlocalized',
                            'scaffold_count_unplaced', 'scaffold_l50', 'scaffold_n50',
                            'contig_count', 'contig_l50', 'contig_n50', 'FtpPath_GenBank']
        self.is_hybrid = True if " x " in self.SpeciesName else False
        self.unknown = True if ("unknown" in self.SpeciesName) or ("Unknown" in self.SpeciesName) else False
        self.only_genus_known = True if (" sp. " in self.SpeciesName) else False
        self.candidate_species = True if ("candidate" in self.SpeciesName) or ("Candidate" in self.SpeciesName) else False


    @staticmethod
    def parse_stats_from_assembly_summary_metadata(metadata):

        end = re.search("/Stats>", metadata).end()
        tmp_dict = xmltodict.parse(metadata[:end])
        stat_dict = OrderedDict()
        for entry in tmp_dict["Stats"]["Stat"]:
            if "scaffold_count" in entry['@category']:
                stat_dict[entry['@category'] + "_" + entry['@sequence_tag']] = int(entry['#text'])
            else:
                stat_dict[entry['@category']] = int(entry['#text'])

        return stat_dict

    def get_header(self):
        return "\t".join(self.header_list)

    def __str__(self, separator="\t"):
        values_list = []
        for key in self.header_list:
            if key == 'WholeGenome':
                values_list.append("Y" if self.full_genome else "N")
            elif key == "ChromosomeLvl":
                values_list.append("Y" if self.chromosome_lvl else "N")
            elif key == "Chloroplast":
                values_list.append("Y" if self.chloroplast else "N")
            elif key == "Mitochondrion":
                values_list.append("Y" if self.chloroplast else "N")
            elif key == "Plasmid":
                values_list.append("Y" if self.plasmid else "N")
            else:
                if isinstance(self[key], str):
                    tmp = (self[key])
                elif isinstance(self[key], list) or isinstance(self[key], set):
                    tmp = ','.join(self[key])
                elif self[key] is None:
                    tmp = "."
                else:
                    tmp = str(self[key])
                values_list.append(tmp)
        #print values_list
        return separator.join(values_list)

    def ambiguous_species(self):
        return True if (self.is_hybrid or self.unknown or self.only_genus_known or self.candidate_species) else False

    def check_integrity(self, min_scaffold_n50=None, min_contig_n50=None, max_scaffold_l50=None,
                        max_contig_l50=None, max_contig_count=None, max_scaffold_count=None,
                        max_chromosome_count=None, min_chromosome_count=None, max_unlocalized_scaffolds=None,
                        max_unplaced_scaffolds=None, max_total_length=None, min_total_length=None,
                        max_ungapped_length=None, min_ungapped_length=None):
        if min_contig_n50:
            if self.contig_n50 < min_contig_n50:
                return False
        if min_scaffold_n50:
            if self.scaffold_n50 < min_scaffold_n50:
                return False
        if max_scaffold_l50:
            if self.scaffold_l50 > max_scaffold_l50:
                return False
        if max_contig_l50:
            if self.contig_l50 > max_contig_l50:
                return False
        if max_contig_count:
            if self.contig_count > max_contig_count:
                return False
        if max_scaffold_count:
            if self.scaffold_count_all > max_scaffold_count:
                return False
        if max_chromosome_count:
            if self.chromosome_count > max_chromosome_count:
                return False
        if min_chromosome_count:
            if self.chromosome_count < min_chromosome_count:
                return False
        if max_unlocalized_scaffolds:
            if self.scaffold_count_unlocalized > max_unlocalized_scaffolds:
                return False
        if max_unplaced_scaffolds:
            if self.scaffold_count_unplaced < max_unplaced_scaffolds:
                return False
        if max_total_length:
            if self.total_length > max_total_length:
                return False
        if min_total_length:
            if self.total_length < min_total_length:
                return False

        if max_ungapped_length:
            if self.ungapped_length > max_ungapped_length:
                return False
        if min_ungapped_length:
            if self.ungapped_length < min_ungapped_length:
                return False

        return True


class AssemblySummaryList(list):
    def __init__(self, entrez_summary_biopython=None):
        list.__init__(self)
        if entrez_summary_biopython:
            for assembly_summary_biopython_dict_element in entrez_summary_biopython['DocumentSummarySet']['DocumentSummary']:
            #print assembly_summary_biopython_dict_element
                self.append(AssemblySummary(assembly_summary_biopython_dict_element))
            if len(self) > 0:
                self.header = self[0].get_header()

        else:
            self.header = "\t".join(['SpeciesName', 'Taxid', 'Isolate', 'Sub_type', 'Sub_value', 'Sex', 'BioSampleId',
                                     'BioSampleAccn', 'WGS', 'AssemblyAccession', 'AssemblyName', 'SubmissionDate',
                                     'LastUpdateDate', 'AssemblyType', 'AssemblyStatus',
                                     'WholeGenome', "ChromosomeLvl", 'Chloroplast', "Mitochondrion", "Plasmid",
                                     'total_length', 'ungapped_length', "alt_loci_count",
                                     'chromosome_count', 'non_chromosome_replicon_count', 'replicon_count',
                                     'scaffold_count_all', 'scaffold_count_placed', 'scaffold_count_unlocalized',
                                     'scaffold_count_unplaced', 'scaffold_l50', 'scaffold_n50',
                                     'contig_count', 'contig_l50', 'contig_n50', 'FtpPath_GenBank'])

    def __str__(self):
        string = self.header + "\n"
        for summary in self:
            string += str(summary) + "\n"

    def write(self, output_file):
        with open(output_file, "w") as out_fd:
            out_fd.write(self.header + "\n")
            for summary in self:
                out_fd.write(str(summary) + "\n")

    def filter(self, expression):
        filtered = AssemblySummaryList()
        filtered_out = AssemblySummaryList()

        for summary in self:
            if expression(summary):
                filtered.append(summary)
            else:
                filtered_out.append(summary)

        return filtered, filtered_out

    def filter_non_chrom_level_genomes(self):
        def chromosome_lvl(summary):
            return summary.chromosome_lvl

        return self.filter(chromosome_lvl)

    def filter_ambiguous_species(self, verbose=False):
        if verbose:
            for summary in self:
                print("\tSpeciesName\tIsAmbigious")
                print("\t%s\t%s" % (str(summary.SpeciesName), str(summary.ambiguous_species())))
        return self.filter(lambda summary: not summary.ambiguous_species())

    def filter_by_integrity(self, min_scaffold_n50=None, min_contig_n50=None, max_scaffold_l50=None,
                            max_contig_l50=None, max_contig_count=None, max_scaffold_count=None,
                            max_chromosome_count=None, min_chromosome_count=None, max_unlocalized_scaffolds=None,
                            max_unplaced_scaffolds=None, max_total_length=None, min_total_length=None,
                            max_ungapped_length=None, min_ungapped_length=None,
                            no_ambiguous_species=False):

        def expression(summary):
            return summary.check_integrity(min_scaffold_n50=min_scaffold_n50, min_contig_n50=min_contig_n50,
                                           max_scaffold_l50=max_scaffold_l50, max_contig_l50=max_contig_l50,
                                           max_contig_count=max_contig_count, max_scaffold_count=max_scaffold_count,
                                           max_chromosome_count=max_chromosome_count,
                                           min_chromosome_count=min_chromosome_count,
                                           max_unlocalized_scaffolds=max_unlocalized_scaffolds,
                                           max_unplaced_scaffolds=max_unplaced_scaffolds,
                                           max_total_length=max_total_length, min_total_length=min_total_length,
                                           max_ungapped_length=max_ungapped_length,
                                           min_ungapped_length=min_ungapped_length) \
                   and ((not summary.ambiguous_species()) if no_ambiguous_species else True)

        return self.filter(expression)


class NCBIRoutines(SequenceRoutines):
    def __init__(self):
        SequenceRoutines.__init__(self)
        self.ncbi_ftp = "ftp://ftp-trace.ncbi.nlm.nih.gov/"

    def get_sra_ftp_path_from_id(self, sra_id):

        sra_reads_dir = "sra/sra-instant/reads/"
        id_code_dict = {
            "DRR": "ByRun",
            "ERR": "ByRun",
            "SRR": "ByRun",

            "DRX": "ByExp",
            "ERX": "ByExp",
            "SRX": "ByExp",

            "DRS": "BySample",
            "ERS": "BySample",
            "SRS": "BySample",

            "DRP": "ByStudy",
            "ERP": "ByStudy",
            "SRP": "ByStudy"
        }
        id_group = sra_id[:3]
        id_subgroup = sra_id[:6]

        id_type = id_code_dict[id_group]

        return "%s%s%s/sra/%s/%s/%s/%s.sra" % (self.ncbi_ftp, sra_reads_dir, id_type, id_group,
                                               id_subgroup, sra_id, sra_id)

    @staticmethod
    def efetch(database, id_list, out_file, retmode=None, rettype=None, seq_start=None, seq_stop=None, strand=None, verbose=False,
               number_of_retries=100, retry_delay=None, log_file="efetch.log"):
        # replacement for Biopython Entrez.efetch
        # Biopython Entrez.efetch is bugged - it ignores seq_start and seq_stop values
        # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=669632474&retmode=text&rettype=gb&seq_start=10832751&seq_stop=10848091&strand=1
        query = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"

        query += "db=%s" % database

        query += "&id=%s" % (id_list if isinstance(id_list, str) else ",".join(id_list))
        query += "&retmode=%s" % retmode if retmode else ""
        query += "&rettype=%s" % rettype if rettype else ""

        query += "&seq_start=%s" % str(seq_start) if seq_start else ""
        query += "&seq_stop=%s" % str(seq_stop) if seq_stop else ""
        query += "&strand=%s" % str(strand) if strand else ""

        curl_options = " --retry %i" % number_of_retries if number_of_retries else ""
        curl_options += " --retry-delay %i" % retry_delay if retry_delay else ""
        curl_options += " '%s'" % query
        curl_options += " > %s" % out_file

        curl_string = "curl %s" % curl_options
        with open(log_file, "a") as log_fd:
            log_fd.write(curl_string + "\n")
        if verbose:
            print(curl_string)
        os.system(curl_string)

    def get_gene_sequences(self, email, query, retmax=100000, output_directory=None):
        if output_directory:
            self.safe_mkdir(output_directory)
        Entrez.email = email
        handle = Entrez.esearch(term=query, db="gene", rettype="xml", retmax=retmax)
        records = Entrez.read(handle, validate=False)
        gene_ids = records['IdList']
        for gene_id in gene_ids:
            print("Extracting sequence for gi %s" % gene_id)
            id_handle = Entrez.efetch("gene", id=gene_id, rettype="xml", retmax=retmax)
            record = Entrez.read(id_handle, validate=False)

            #print(record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"])
            if "Entrezgene_locus" not in record[0]:
                continue
            if "Gene-commentary_seqs" not in record[0]["Entrezgene_locus"][0]:
                continue
            if "Seq-loc_int" not in record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]:
                continue
            gi = record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_id']['Seq-id']['Seq-id_gi']
            start = int(record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_from'])
            end = int(record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_to']) + 1
            strand = record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_strand']['Na-strand'].attributes['value']
            strand = 1 if strand == "plus" else 2
            # print(gi)
            # print(start, end, strand)
            # Biopython Entrez.efetch is bugged - it ignores seq_start and seq_stop values
            # eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=669632474&retmode=text&rettype=gb&seq_start=10832751&seq_stop=10848091&strand=1

            out_file = "%s/%s.gb" % (output_directory, gi)

            self.efetch("nuccore", gi, out_file, retmode="text", rettype="gb", seq_start=start, seq_stop=end,
                        strand=strand)
            time.sleep(0.4)

    @staticmethod
    def get_longest_proteins_from_protein_tab_file(protein_tab_file, output_prefix):
        sorted_tab_file = "%s.sorted.tab" % output_prefix
        longest_pep_tab_file = "%s.longest_pep.tab" % output_prefix
        longest_pep_id_file = "%s.longest_pep.ids" % output_prefix

        with open(protein_tab_file)as tmp_fd:
            tmp = tmp_fd.readline().strip().split("\t")
            gene_id_index = tmp.index("GeneID")
            protein_id_index = tmp.index("Protein product")
            length_index = tmp.index("Length")

        sort_string = "(head -n 1 %s && tail -n +2 %s | sort -k%i,%i -k%i,%inr) > %s" % (protein_tab_file,
                                                                                         protein_tab_file,
                                                                                         gene_id_index + 1,
                                                                                         gene_id_index + 1,
                                                                                         length_index + 1,
                                                                                         length_index + 1,
                                                                                         sorted_tab_file)

        os.system(sort_string)

        longest_pep_ids = IdSet()

        with open(sorted_tab_file, "r") as in_fd:
            with open(longest_pep_tab_file, "w") as out_fd:
                out_fd.write(in_fd.readline())
                prev_gene = ""
                for line in in_fd:
                    tmp = line.split("\t")
                    if tmp[gene_id_index] == prev_gene:
                        continue
                    prev_gene = tmp[gene_id_index]
                    out_fd.write(line)
                    longest_pep_ids.add(tmp[protein_id_index])

        longest_pep_ids.write(longest_pep_id_file)

    def get_cds_for_proteins(self, protein_id_list, output_prefix, download_chunk_size=100, temp_dir_prefix="temp",
                             keyword_for_mitochondrial_pep="mitochondrion"):
        duplicated_protein_id_list = IdList()
        protein_id_uniq_list = IdList()

        prev_protein_id = ""
        for protein_id in sorted(protein_id_list):
            if protein_id != prev_protein_id:
                protein_id_uniq_list.append(protein_id)
                prev_protein_id = protein_id

            else:
                duplicated_protein_id_list.append(protein_id)

        duplicated_protein_id_list.write("%s.duplicated.pep.ids" % output_prefix)
        if duplicated_protein_id_list:
            print("WARNING!!! Duplicated protein ids were detected...")
            with self.metaopen("%s.warnings" % output_prefix, "w") as warn_fd:
                warn_fd.write("WARNING!!! Duplicated protein ids were detected...")

        from RouToolPa.Tools.Abstract import Tool
        output_directory = self.check_dir_path(self.split_filename(output_prefix)[0])
        transcript_temp_dir = "%s%s_transcripts" % (output_directory, temp_dir_prefix)
        protein_temp_dir = "%s%s_proteins" % (output_directory, temp_dir_prefix)
        number_of_ids = len(protein_id_uniq_list)

        print("Total %i ids" % number_of_ids)

        for directory in transcript_temp_dir, protein_temp_dir:
            print(directory)
            self.safe_mkdir(directory)
        pep_file = "%s.pep.genbank" % output_prefix
        transcript_file = "%s.trascript.genbank" % output_prefix
        mito_transcript_file = "%s.trascript.mito.genbank" % output_prefix

        ranges = np.append(np.arange(0, number_of_ids, download_chunk_size), [number_of_ids])

        print("Downloading proteins...")
        for i in range(0, len(ranges)-1):
            print("Downloading chunk %i" % i)
            pep_tmp_file = "%s/chunk_%i" % (protein_temp_dir, i)
            self.efetch("protein", protein_id_uniq_list[ranges[i]:ranges[i+1]], pep_tmp_file, rettype="gb", retmode="text",
                        log_file="%s.pep.efetch.log" % output_prefix)

        os.system("cat %s/* > %s" % (protein_temp_dir, pep_file))

        peptide_dict = SeqIO.index_db("%stmp.idx" % output_directory, pep_file, format="genbank")
        downloaded_protein_ids = IdList(peptide_dict.keys())

        print("%i proteins were downloaded" % len(downloaded_protein_ids))
        not_downloaded_proteins_ids = Tool.intersect_ids([protein_id_uniq_list], [downloaded_protein_ids], mode="only_a")
        print("%i proteins were not downloaded" % len(not_downloaded_proteins_ids))
        not_downloaded_proteins_ids.write("%s.not_downloaded.ids" % output_prefix)
        downloaded_protein_ids.write("%s.downloaded.ids" % output_prefix)
        print(Tool.intersect_ids([protein_id_uniq_list], [downloaded_protein_ids], mode="count"))

        pep_without_transcripts = IdList()
        mito_pep_without_transcripts = IdList()
        pep_with_several_CDS_features = IdList()
        mito_pep_with_several_CDS_features = IdList()
        pep_to_transcript_accordance = SynDict()
        mito_pep_to_transcript_accordance = SynDict()
        transcript_ids = IdList()
        mito_transcript_ids = IdList()
        print("Extracting transcript ids corresponding to proteins...")

        for pep_id in peptide_dict:
            if keyword_for_mitochondrial_pep in peptide_dict[pep_id].annotations["source"]:
                for feature in peptide_dict[pep_id].features:
                    if feature.type == "CDS":
                        try:
                            transcript_id = feature.qualifiers["coded_by"][0].split(":")[0]
                            if pep_id not in mito_pep_to_transcript_accordance:
                                mito_pep_to_transcript_accordance[pep_id] = [transcript_id]
                            else:
                                mito_pep_to_transcript_accordance[pep_id].append(transcript_id)
                                print("Genbank record for %s contains several CDS features" % pep_id)
                                mito_pep_with_several_CDS_features.append(pep_id)
                            if transcript_id in mito_transcript_ids:
                                print("Repeated transcript id: %s" % transcript_id)
                                continue
                            mito_transcript_ids.append(transcript_id)
                        except:
                            print("Transcript id for %s was not found" % pep_id)
                            mito_pep_without_transcripts.append(pep_id)

            else:
                for feature in peptide_dict[pep_id].features:
                    if feature.type == "CDS":
                        try:
                            transcript_id = feature.qualifiers["coded_by"][0].split(":")[0]
                            if pep_id not in pep_to_transcript_accordance:
                                pep_to_transcript_accordance[pep_id] = [transcript_id]
                            else:
                                pep_to_transcript_accordance[pep_id].append(transcript_id)
                                print("Genbank record for %s contains several CDS features" % pep_id)
                                pep_with_several_CDS_features.append(pep_id)
                            if transcript_id in transcript_ids:
                                print("Repeated transcript id: %s" % transcript_id)
                                continue
                            transcript_ids.append(transcript_id)
                        except:
                            print("Transcript id for %s was not found" % pep_id)
                            pep_without_transcripts.append(pep_id)

        pep_with_several_CDS_features.write("%s.pep_with_several_CDS.ids" % output_prefix)
        pep_without_transcripts.write("%s.pep_without_transcripts.ids" % output_prefix)
        transcript_ids.write("%s.transcripts.ids" % output_prefix)

        mito_pep_with_several_CDS_features.write("%s.mito.pep_with_several_CDS.ids" % output_prefix)
        mito_pep_without_transcripts.write("%s.mito.pep_without_transcripts.ids" % output_prefix)
        mito_transcript_ids.write("%s.mito.transcripts.ids" % output_prefix)

        number_of_transcripts = len(transcript_ids) + len(mito_transcript_ids)
        print("%i transcripts were found" % number_of_transcripts)

        pep_to_transcript_accordance.write("%s.pep_to_transcript.accordance" % output_prefix, splited_values=True)
        mito_pep_to_transcript_accordance.write("%s.mito.pep_to_transcript.accordance" % output_prefix, splited_values=True)

        transcript_ranges = np.append(np.arange(0, number_of_transcripts, download_chunk_size), [number_of_transcripts])

        print("Downloading transcripts...")
        for i in range(0, len(transcript_ranges)-1):
            print("Downloading chunk %i" % i)
            transcript_tmp_file = "%s/chunk_%i" % (transcript_temp_dir, i)
            self.efetch("nuccore", transcript_ids[transcript_ranges[i]:transcript_ranges[i+1]],
                        transcript_tmp_file, rettype="gb", retmode="text",
                        log_file="%s.transcripts.efetch.log" % output_prefix)

        os.system("cat %s/* > %s" % (transcript_temp_dir, transcript_file))

        transcript_dict = SeqIO.index_db("%stmp_1.idx" % output_directory, transcript_file, format="genbank")

        print("Downloading mitochondrial transcripts...")
        self.efetch("nuccore", mito_transcript_ids, mito_transcript_file, rettype="gb", retmode="text",
                    log_file="%s.mito.transcripts.efetch.log" % output_prefix)
        mito_transcript_dict = SeqIO.index_db("%stmp_2.idx" % output_directory, mito_transcript_file, format="genbank")

        actual_protein_ids_from_transcripts = IdList()
        actual_mito_protein_ids_from_transcripts = IdList()
        actual_absent_protein_ids_from_transcripts = IdList()
        actual_absent_mito_protein_ids_from_transcripts = IdList()
        actual_pep_to_transcript_accordance = SynDict()
        actual_mito_pep_to_transcript_accordance = SynDict()
        actual_transcript_ids = IdList()
        actual_mito_transcript_ids = IdList()

        cds_records_dict = OrderedDict()
        for transcript_id in transcript_dict:
            CDS_number = 0
            for feature in transcript_dict[transcript_id].features:
                if feature.type == "CDS":
                    CDS_number += 1
            if CDS_number > 1:
                CDS_counter = 1
                for feature in transcript_dict[transcript_id].features:
                    if feature.type == "CDS":
                        feature_seq = feature.extract(transcript_dict[transcript_id].seq)
                        feature_id = "%s.%i" % (transcript_id, CDS_counter)
                        if "protein_id" in feature.qualifiers:
                            description = "protein=%s" % feature.qualifiers["protein_id"][0]
                            actual_protein_ids_from_transcripts.append(feature.qualifiers["protein_id"][0])
                            actual_pep_to_transcript_accordance[feature.qualifiers["protein_id"][0]] = [feature_id]
                            actual_transcript_ids.append(feature_id)
                        else:
                            print("Corresponding protein id was not found for %s " % transcript_id)

                        cds_records_dict[feature_id] = SeqRecord(seq=feature_seq, id=feature_id, description=description)
                        CDS_counter += 1
            else:
                for feature in transcript_dict[transcript_id].features:
                    if feature.type == "CDS":
                        feature_seq = feature.extract(transcript_dict[transcript_id].seq)
                        feature_id = transcript_id
                        if "protein_id" in feature.qualifiers:
                            description = "protein=%s" % feature.qualifiers["protein_id"][0]
                            actual_protein_ids_from_transcripts.append(feature.qualifiers["protein_id"][0])
                            actual_pep_to_transcript_accordance[feature.qualifiers["protein_id"][0]] = [feature_id]
                            actual_transcript_ids.append(feature_id)
                        else:
                            print("Corresponding protein id was not found for %s " % transcript_id)

                        cds_records_dict[feature_id] = SeqRecord(seq=feature_seq, id=feature_id, description=description)
                        break
        SeqIO.write(self.record_by_id_generator(cds_records_dict), "%s.all.cds" % output_prefix, format="fasta")
        SeqIO.write(self.record_by_id_generator(peptide_dict, id_list=actual_protein_ids_from_transcripts),
                    "%s.with_cds.pep" % output_prefix, format="fasta")
        SeqIO.write(self.record_by_id_generator(cds_records_dict, id_list=actual_transcript_ids),
                    "%s.with_pep.cds" % output_prefix, format="fasta")

        actual_pep_to_transcript_accordance.write("%s.pep_to_transcript.actual.accordance" % output_prefix, splited_values=True)

        mito_cds_records_dict = OrderedDict()
        for transcript_id in mito_transcript_dict:
            CDS_number = 0
            #print transcript_id
            for feature in mito_transcript_dict[transcript_id].features:
                if feature.type == "CDS":
                    CDS_number += 1
            #print CDS_number
            if CDS_number > 1:
                CDS_counter = 1
                for feature in mito_transcript_dict[transcript_id].features:
                    if feature.type == "CDS":
                        feature_seq = feature.extract(mito_transcript_dict[transcript_id].seq)
                        feature_id = "%s.%i" % (transcript_id, CDS_counter)
                        if "protein_id" in feature.qualifiers:
                            description = "protein=%s" % feature.qualifiers["protein_id"][0]
                            actual_mito_protein_ids_from_transcripts.append(feature.qualifiers["protein_id"][0])
                            actual_mito_pep_to_transcript_accordance[feature.qualifiers["protein_id"][0]] = [feature_id]
                            actual_mito_transcript_ids.append(feature_id)
                        else:
                            print("Corresponding protein id was not found for %s " % transcript_id)
                        CDS_counter += 1
                        mito_cds_records_dict[feature_id] = SeqRecord(seq=feature_seq, id=feature_id, description=description)

            else:
                for feature in mito_transcript_dict[transcript_id].features:
                    if feature.type == "CDS":
                        feature_seq = feature.extract(mito_transcript_dict[transcript_id].seq)
                        feature_id = transcript_id
                        if "protein_id" in feature.qualifiers:
                            description = "protein=%s" % feature.qualifiers["protein_id"][0]
                            actual_mito_protein_ids_from_transcripts.append(feature.qualifiers["protein_id"][0])
                            actual_mito_pep_to_transcript_accordance[feature.qualifiers["protein_id"][0]] = [feature_id]
                            actual_mito_transcript_ids.append(feature_id)
                        else:
                            print("Corresponding protein id was not found for %s " % transcript_id)

                        mito_cds_records_dict[feature_id] = SeqRecord(seq=feature_seq, id=feature_id, description=description)
                        break
        SeqIO.write(self.record_by_id_generator(mito_cds_records_dict), "%s.mito.all.cds" % output_prefix, format="fasta")
        SeqIO.write(self.record_by_id_generator(peptide_dict, id_list=actual_mito_protein_ids_from_transcripts),
                    "%s.mito.with_cds.pep" % output_prefix, format="fasta")
        SeqIO.write(self.record_by_id_generator(mito_cds_records_dict, id_list=actual_mito_transcript_ids),
                    "%s.mito.with_pep.cds" % output_prefix, format="fasta")

        #print actual_mito_pep_to_transcript_accordance
        actual_mito_pep_to_transcript_accordance.write("%s.mito.pep_to_transcript.actual.accordance" % output_prefix, splited_values=True)

        stat_string = "Input protein ids\t %i\n" % number_of_ids
        stat_string += "Downloaded proteins\t%i\n" % len(downloaded_protein_ids)
        stat_string += "Downloaded transcripts\t%i\n" % (len(transcript_dict) + len(mito_transcript_dict))
        stat_string += "Valid protein and CDS pair\t%i\n" % (len(actual_pep_to_transcript_accordance) +
                                                             len(actual_mito_pep_to_transcript_accordance))
        print(stat_string)

        with open("%s.stats" % output_prefix, "w") as stat_fd:
            stat_fd.write(stat_string)

        description_string = """DESCRIPTION OF FILES CREATED DURING DOWNLOADa AND EXTRACTION OF CDS FOR PROTEIN IDS FROM NCBI

* - prefix of files

RESULTING FILES
*.with_cds.pep                                 Final protein set for nuclear genes
*.with_pep.cds                                 Final CDS set for nuclear genes
*.pep_to_transcript.actual.accordance          Correspondence file between proteins and CDSs for nuclear genes

*.mito.with_cds.pep                            Final protein set for mitochondrial genes
*.mito.with_pep.cds                            Final CDS set for mitochondrial genes
*.mito.pep_to_transcript.actual.accordance     Correspondence file between proteins
                                               and CDSs for mitochondrial genes

INTERMIDIATE FILES AND DIRECTORIES

temp_proteins/                                 Directory with downloaded proteins
temp_transcripts/                              Directory with downloaded transcripts

*.pep.genbank                                  Merged file with downloaded nucleolar proteins

*.downloaded.ids                               Ids of downloaded proteins
*.not_downloaded.ids                           Ids of not downloaded proteins

*.trascript.genbank                            Merged file with downloaded nucleolar transcripts
*.trascript.mito.genbank                       Merged file with downloaded mitochondrial
                                               transcripts

*.mito.all.cds
*.mito.pep_to_transcript.accordance
*.mito.pep_without_transcripts.ids
*.mito.pep_with_several_CDS.ids
*.mito.transcripts.ids

*.all.cds
*.pep_to_transcript.accordance
*.pep_without_transcripts.ids
*.pep_with_several_CDS.ids

*.stats                                       Statistics
*.transcripts.ids
"""

        with open("%sDESCRIPTION" % output_directory, "w") as descr_fd:
            descr_fd.write(description_string)

        for filename in "tmp.idx", "tmp_1.idx", "tmp_2.idx":
            os.remove("%s%s" % (output_directory, filename))

    def get_cds_for_proteins_from_id_file(self, protein_id_file, output_prefix, download_chunk_size=100):
        pep_ids = IdList()
        pep_ids.read(protein_id_file)

        self.get_cds_for_proteins(pep_ids, output_prefix, download_chunk_size=download_chunk_size)

    @staticmethod
    def safe_entrez_function(entrez_function, max_tries=10, **kwargs):

        for entry in range(0, max_tries):
            try:
                result = entrez_function(**kwargs) #Entrez.read(Entrez.efetch(db="taxonomy", id=taxon, retmode="xml"))
                if result:
                    return result
            except:
                pass

        return None

    def get_taxonomy(self, taxa_list, output_file, email, input_type="latin", max_tries=10):
        Entrez.email = email
        out_file = open(output_file, "w")
        out_file.write("#species\trank\tlatin_name\tlineage\n")

        species_syn_dict = SynDict()

        if input_type == "latin":
            for taxon in taxa_list:
                print("Handling %s" % taxon)
                summary = Entrez.read(Entrez.esearch(db="taxonomy", term=taxon))
                if summary:
                    id_list = summary["IdList"]
                    species_syn_dict[taxon] = []
                    for id in id_list:
                        print("\tHandling %s" % id)
                        record = self.safe_entrez_function(Entrez.efetch, max_tries=max_tries, db="taxonomy", id=id, retmode="xml")
                        if record:
                            record = Entrez.read(record)
                        else:
                            continue
                        
                        out_file.write("%s\t%s\t%s\t%s\n" % (taxon,
                                                             record[0]["Rank"],
                                                             record[0]['ScientificName'],
                                                             record[0]["Lineage"]))

                        species_syn_dict[taxon].append(record[0]['ScientificName'])
                        #species_set.add(record[0]["Species"])
        elif input_type == "id":
            for taxon in taxa_list:
                print("Handling %s" % taxon)
                species_syn_dict[taxon] = []
                #print taxon

                record = self.safe_entrez_function(Entrez.efetch, db="taxonomy", id=taxon, retmode="xml")
                if record:
                    record = Entrez.read(record)
                else:
                    continue
                #print record
                out_file.write("%s\t%s\t%s\t%s\n" % (taxon,
                                                     record[0]["Rank"],
                                                     record[0]['ScientificName'],
                                                     record[0]["Lineage"]))
                #print record[0]
                species_syn_dict[taxon].append(record[0]['ScientificName'])
                #species_set.add(record[0]["Species"])
        out_file.close()
        return species_syn_dict

    @staticmethod
    def get_subtaxa_for_taxa(taxa_list, email, output_prefix, subtaxa_rank=('species'), verbose=False):
        Entrez.email = email
        taxonomy_file = "%s.taxonomy" % output_prefix
        species_syn_file = "%s.species.ids" % output_prefix

        subtaxa_rank_list = subtaxa_rank.split(",") if isinstance(subtaxa_rank, str) else subtaxa_rank
        species_syn_dict = SynDict()
        with open(taxonomy_file, "w") as taxa_fd:
            taxa_fd.write("#species\trank\tlatin_name\tcommon_name\tlineage\n")
            for taxon in taxa_list:
                print("Handling %s..." % taxon)
                summary = Entrez.read(Entrez.esearch(db="taxonomy",
                                                     term="%s[subtree] AND (%s)" % (taxon, "[rank] OR ".join(subtaxa_rank_list) if len(subtaxa_rank_list) > 4 else "%s[rank]" % subtaxa_rank_list[0]), retmax=4000))

                if verbose:
                    print (summary)

                if summary:
                    id_list = summary["IdList"]
                    for id in id_list:
                        print("handling %s..." % id)
                        record = Entrez.read(Entrez.efetch(db="taxonomy", id=id, retmode="xml"))
                        if verbose:
                            print(record)
                        common_name = "NA"
                        if 'OtherNames' in record[0]:
                            if "GenbankCommonName" in record[0]['OtherNames']:
                                common_name = record[0]['OtherNames']["GenbankCommonName"]

                        taxa_fd.write("%s\t%s\t%s\t%s\n" % (record[0]['ScientificName'],
                                                                record[0]["Rank"],
                                                                common_name,
                                                                record[0]["Lineage"]))

                        species_syn_dict[record[0]['ScientificName']] = common_name

        species_syn_dict.write(filename=species_syn_file, header="#latin_name\tcommon_name\n")
        return species_syn_dict

    def get_taxonomy_from_id_file(self, taxa_file, output_file, email, input_type="latin", column=0, separator="\t"):

        taxa_list = IdList(filename=taxa_file, column_number=column, column_separator=separator)

        return self.get_taxonomy(taxa_list, output_file, email, input_type=input_type)

    def get_protein_info_from_ncbi_gff(self, ncbi_gff, output_prefix):
        ncbi_gff_collection = CollectionGFF
        for record in ncbi_gff_collection.gff_simple_generator(ncbi_gff):
            pass

    @staticmethod
    def extract_gene_info_from_ncbi_gff(input_gff, output, feature_type_list=["gene"]):
        with open(input_gff, "r") as in_fd:
            with open(output, 'w') as out_fd:
                out_fd.write("#GeneID\tSource\tName\tGeneSynonym\tGeneBiotype\tDbxref\tDescription\n")
                for line in in_fd:
                    if line[0] == "#":
                        continue
                    tmp_list = line.strip().split("\t")
                    if tmp_list[2] not in feature_type_list:
                        continue
                    source = tmp_list[1]

                    annotation_dict = AnnotationsRoutines.parse_gff_annotation_string_to_dict(tmp_list[8])
                    if "Dbxref" in annotation_dict:
                        for entry in annotation_dict["Dbxref"]:
                            if "GeneID" in entry:
                                geneid = entry.split(":")[1]
                                break
                            else:
                                geneid = "."
                        dbxref = ",".join(annotation_dict["Dbxref"])
                    else:
                        geneid = "."
                        dbxref = "."
                    name = ",".join(annotation_dict["Name"] if "Name" in annotation_dict else ".")
                    description = ",".join(annotation_dict["description"] if "description" in annotation_dict else ".")
                    gene_biotype = ",".join(annotation_dict["gene_biotype"] if "gene_biotype" in annotation_dict else ".")
                    gene_synonym = ",".join(annotation_dict["gene_synonym"] if "gene_synonym" in annotation_dict else ".")
                    out_fd.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneid, source, name, gene_synonym,
                                                                   gene_biotype, dbxref, description))

    def get_taxa_genomes_summary(self, taxa, email, output_directory, output_prefix,
                                 max_ids_per_query=8000, max_download_attempts=500,
                                 min_scaffold_n50=None, min_contig_n50=None, max_scaffold_l50=None,
                                 max_contig_l50=None, max_contig_count=None, max_scaffold_count=None,
                                 max_chromosome_count=None, min_chromosome_count=None, max_unlocalized_scaffolds=None,
                                 max_unplaced_scaffolds=None, max_total_length=None, min_total_length=None,
                                 max_ungapped_length=None, min_ungapped_length=None,
                                 no_ambiguous_species=True):
        self.safe_mkdir(output_directory)
        Entrez.email = email
        taxa_list = taxa if isinstance(taxa, Iterable) else [taxa]

        all_files_dir = "%s%s/" % (self.check_path(output_directory), "all")
        nonambiguous_species_all_dir = "%snonambiguous_species_all/" % self.check_path(output_directory)
        ambiguous_species_all_dir = "%s%s/" % (self.check_path(output_directory), "ambiguous_species_all")
        chromosome_lvl_dir = "%s%s/" % (self.check_path(output_directory), "chromosome_lvl")
        non_chromosome_lvl_dir = "%s%s/" % (self.check_path(output_directory), "nonchromosome_lvl")

        filtered_by_integrity_dir = "%s%s/" % (self.check_path(output_directory), "passed_integrity_filters")
        filtered_out_by_integrity_dir = "%s%s/" % (self.check_path(output_directory), "not_passed_integrity_filters")

        stat_dir = "%s%s/" % (self.check_path(output_directory), "stat")
        taxa_stat_dir = "%s%s/" % (self.check_path(output_directory), "taxa_stat")
        for subdir in (all_files_dir, chromosome_lvl_dir, non_chromosome_lvl_dir, stat_dir,
                       taxa_stat_dir, nonambiguous_species_all_dir, ambiguous_species_all_dir):
            self.safe_mkdir(subdir)

        filter_by_integrity = min_scaffold_n50 or min_contig_n50 or max_scaffold_l50 or max_contig_l50 \
                              or max_contig_count or max_scaffold_count or max_chromosome_count \
                              or min_chromosome_count or max_unlocalized_scaffolds \
                              or max_unplaced_scaffolds or max_total_length or min_total_length \
                              or max_ungapped_length or min_ungapped_length

        if filter_by_integrity:
            for subdir in (filtered_by_integrity_dir, filtered_out_by_integrity_dir):
                self.safe_mkdir(subdir)

        for taxon in taxa_list:
            search_term = "%s[Orgn]" % taxon

            attempt_counter = 1
            while True:
                try:
                    summary = Entrez.read(Entrez.esearch(db="genome", term=search_term, retmax=10000, retmode="xml"))
                    break
                except URLError:
                    if attempt_counter > max_download_attempts:
                        URLError("Network problems. Maximum attempt number is exceeded")
                    print("URLError. Retrying... Attempt %i" % attempt_counter)
                    attempt_counter += 1

            print("Were found %s species" % summary["Count"])
            #print summary

            taxon_stat_file = "%s/%s.stat" % (taxa_stat_dir, taxon.replace(" ", "_"))
            taxon_stat_dict = TwoLvlDict()

            for species_id in summary["IdList"]: #[167] :
                print("Handling species id %s " % species_id)

                species_stat_file = "%s/%s.stat" % (stat_dir, species_id)
                species_stat_dict = TwoLvlDict()
                species_stat_dict[species_id] = OrderedDict()

                taxon_stat_dict[species_id] = OrderedDict()

                for stat in "all", "chromosome_lvl", "non_chromosome_lvl":
                    species_stat_dict[species_id][stat] = 0
                    taxon_stat_dict[species_id][stat] = 0
                #species_summary = Entrez.read(Entrez.esummary(db="genome", id=species_id, retmax=10000, retmode="xml"))
                #print species_summary

                # get assemblies linked with genome of species

                attempt_counter = 1
                while True:
                    try:
                        assembly_links = Entrez.read(Entrez.elink(dbfrom="genome", id=species_id, retmode="xml",
                                                                  retmax=10000, linkname="genome_assembly"))
                        break
                    except URLError:
                        if attempt_counter > max_download_attempts:
                            URLError("Network problems. Maximum attempt number is exceeded")
                        print("URLError. Retrying... Attempt %i" % attempt_counter)
                        attempt_counter += 1

                assembly_number = len(assembly_links)
                #print links
                #print links[0]["LinkSetDb"][0]["Link"]
                if assembly_links:
                    if "LinkSetDb" in assembly_links[0]:
                        if assembly_links[0]["LinkSetDb"]:
                            if "Link" in assembly_links[0]["LinkSetDb"][0]:
                                assembly_ids = [id_dict["Id"] for id_dict in assembly_links[0]["LinkSetDb"][0]["Link"]]
                            else:
                                continue
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
                number_of_ids = len(assembly_ids)

                print("\tFound %i assemblies" % number_of_ids)

                id_group_edges = np.arange(0, number_of_ids+1, max_ids_per_query)

                if id_group_edges[-1] != number_of_ids:
                    id_group_edges = np.append(id_group_edges, number_of_ids)

                number_of_id_groups = len(id_group_edges) - 1

                #print len(assembly_links[0]["LinkSetDb"][0]["Link"])
                #print assembly_ids
                #print len(assembly_ids)
                #assembly_dict = TwoLvlDict()
                #assemblies_with_ambiguous_taxonomies = SynDict()
                #summaries = Entrez.read(Entrez.esummary(db="assembly", id=",".join(assembly_ids), retmode="xml"))

                summary_list = None
                for i in range(0, number_of_id_groups):
                    print("\tDownloading summary about assemblies %i - %i" % (id_group_edges[i]+1, id_group_edges[i+1]))
                    #print len(assembly_ids[id_group_edges[i]:id_group_edges[i+1]])
                    summaries = Entrez.read(Entrez.esummary(db="assembly",
                                                            id=",".join(assembly_ids[id_group_edges[i]:id_group_edges[i+1]]),
                                                            retmode="xml"), validate=False)
                    tmp_summary_list = AssemblySummaryList(entrez_summary_biopython=summaries)
                    summary_list = (summary_list + tmp_summary_list) if summary_list else tmp_summary_list

                print("\tDownloaded %i" % len(summary_list))

                if len(summary_list) != number_of_ids:
                    print("\tWARNING:Not all assemblies were downloaded")
                    """
                    print("\tFollowing assemblies were not downloaded(ids):%s" % ",".join(set()))
                    """

                if summary_list:
                    species_stat_dict[species_id]["all"] = len(summary_list)
                    taxon_stat_dict[species_id]["all"] = len(summary_list)
                    output_file = "%s%s.genome.summary" % ((output_prefix + ".") if output_prefix else "", species_id)
                                                           #summary_list[0]['SpeciesName'].replace(" ", "_"))

                    all_output_file = "%s/%s" % (all_files_dir, output_file)
                    chromosome_lvl_output_file = "%s/%s" % (chromosome_lvl_dir, output_file)
                    non_chromosome_lvl_output_file = "%s/%s" % (non_chromosome_lvl_dir, output_file)
                    nonambiguous_species_output_file = "%s/%s" % (nonambiguous_species_all_dir, output_file)
                    ambiguous_species_output_file = "%s/%s" % (ambiguous_species_all_dir, output_file)
                    chromosome_lvl_summary_list, non_chromosome_lvl_summary_list = summary_list.filter_non_chrom_level_genomes()
                    filtered_by_integrity_file = "%s/%s" % (filtered_by_integrity_dir, output_file)
                    filtered_out_by_integrity_file = "%s/%s" % (filtered_out_by_integrity_dir, output_file)

                    species_stat_dict[species_id]["chromosome_lvl"] = len(chromosome_lvl_summary_list)
                    taxon_stat_dict[species_id]["chromosome_lvl"] = len(chromosome_lvl_summary_list)
                    species_stat_dict[species_id]["non_chromosome_lvl"] = len(non_chromosome_lvl_summary_list)
                    taxon_stat_dict[species_id]["non_chromosome_lvl"] = len(non_chromosome_lvl_summary_list)

                    print("\tChromosome level assemblies %i" % species_stat_dict[species_id]["chromosome_lvl"])
                    print("\tNon chromosome level assemblies %i" % species_stat_dict[species_id]["non_chromosome_lvl"])

                    if chromosome_lvl_summary_list:
                        chromosome_lvl_summary_list.write(chromosome_lvl_output_file)

                    if non_chromosome_lvl_summary_list:
                        non_chromosome_lvl_summary_list.write(non_chromosome_lvl_output_file)

                    nonambiguous_species_summary_list, ambiguous_species_summary_list = summary_list.filter_ambiguous_species()
                    #print(len(nonambiguous_species_summary_list), len(ambiguous_species_summary_list))
                    species_stat_dict[species_id]["nonambiguous_species"] = len(nonambiguous_species_summary_list)
                    species_stat_dict[species_id]["ambiguous_species"] = len(ambiguous_species_summary_list)
                    print("\tAmbiguous species %i" % species_stat_dict[species_id]["ambiguous_species"])
                    if nonambiguous_species_summary_list:
                        nonambiguous_species_summary_list.write(nonambiguous_species_output_file)
                    if ambiguous_species_summary_list:
                        ambiguous_species_summary_list.write(ambiguous_species_output_file)

                    summary_list.write(all_output_file)

                    if filter_by_integrity:
                        filtered_by_integrity, filtered_out_by_integrity = summary_list.filter_by_integrity(min_scaffold_n50=min_scaffold_n50,
                                                                                                            min_contig_n50=min_contig_n50,
                                                                                                            max_scaffold_l50=max_scaffold_l50,
                                                                                                            max_contig_l50=max_contig_l50,
                                                                                                            max_contig_count=max_contig_count,
                                                                                                            max_scaffold_count=max_scaffold_count,
                                                                                                            max_chromosome_count=max_chromosome_count,
                                                                                                            min_chromosome_count=min_chromosome_count,
                                                                                                            max_unlocalized_scaffolds=max_unlocalized_scaffolds,
                                                                                                            max_unplaced_scaffolds=max_unplaced_scaffolds,
                                                                                                            max_total_length=max_total_length,
                                                                                                            min_total_length=min_total_length,
                                                                                                            max_ungapped_length=max_ungapped_length,
                                                                                                            min_ungapped_length=min_ungapped_length,
                                                                                                            no_ambiguous_species=no_ambiguous_species)
                        species_stat_dict[species_id]["filtered_by_integrity"] = len(filtered_by_integrity)
                        species_stat_dict[species_id]["filtered_out_by_integrity"] = len(filtered_out_by_integrity)
                        if filtered_by_integrity:
                            filtered_by_integrity.write(filtered_by_integrity_file)
                        if filtered_out_by_integrity:
                            filtered_out_by_integrity.write(filtered_out_by_integrity_file)
                        print("\tPassed integrity filters %i" % species_stat_dict[species_id]["filtered_by_integrity"])
                species_stat_dict.write(species_stat_file)

                print("\n\n")

            taxon_stat_dict.write(taxon_stat_file)

            """
                for tmp in summary_set:
                    #print "aaaaaaaaaaa"

                    tmp = Entrez.read(Entrez.esummary(db="assembly",
                                                                                   id=assembly_id,
                                                                            retmode="xml"))

                    #print tmp
                    #assembly_summary = AssemblySummary(tmp)
                    #print assembly_summary

                    #print assembly_summary.SpeciesName


                    taxonomy_links = Entrez.read(Entrez.elink(dbfrom="assembly", id=assembly_id, retmode="xml",
                                                              retmax=10000, linkname="assembly_taxonomy"))
                    #print taxonomy_links
                    taxonomy_ids = [id_dict["Id"] for id_dict in taxonomy_links[0]["LinkSetDb"][0]["Link"]]

                    if len(taxonomy_ids) > 1:
                        print( "WARNING: more than one taxonomy id for assembly %s" % assembly_ids)
                        assemblies_with_ambiguous_taxonomies[assembly_id] = taxonomy_ids
                    #print taxonomy_ids

                    taxonomy_record = Entrez.read(Entrez.efetch(db="taxonomy", id=taxonomy_ids[0], retmode="xml"))

                    #print taxonomy_record
                    scientific_name = taxonomy_record[0]["ScientificName"]
                    lineage = taxonomy_record[0]["Lineage"]
                    #common_name_list = taxonomy_record[0]["CommonName"]
                    print (taxonomy_record[0]["ScientificName"])
                """


"""
assembly_nuccore 	Nucleotide 	Nucleotide 	Nucleotide 	5000
assembly_nuccore_insdc 	Nucleotide INSDC 	Nucleotide INSDC 	Nucleotide INSDC (GenBank) 	5000
assembly_nuccore_refseq 	Nucleotide RefSeq 	Nucleotide RefSeq 	Nucleotide RefSeq 	5000
assembly_nuccore_wgsmaster 	WGS Master 	WGS Master 	WGS Master 	5000
assembly_nuccore_wgscontig 	WGS contigs 	WGS contigs 	WGS contigs 	5000
assembly_nuccore_insdc 	Nucleotide INSDC 	Nucleotide INSDC 	Nucleotide INSDC (GenBank) 	5000
assembly_nuccore_refseq 	Nucleotide RefSeq 	Nucleotide RefSeq 	Nucleotide RefSeq 	5000
assembly_nuccore_wgsmaster 	WGS Master 	WGS Master 	WGS Master
"""

"""

#Example of assembly summary

{u'DocumentSummarySet':
     DictElement({u'DbBuild': 'Build161130-1600.1',
                  u'DocumentSummary': [DictElement({u'ChainId': '1886755',
                                                    u'AsmUpdateDate': '2016/11/28 00:00',
                                                    u'GbUid': '3745278',
                                                    u'PropertyList': ['full-genome-representation', 'has-chromosome', 'has-plasmid', 'latest', 'latest_genbank'],
                                                    u'SubmissionDate': '2016/11/28 00:00',
                                                    u'GB_BioProjects': [{u'BioprojectAccn': 'PRJNA353939',
                                                                         u'BioprojectId': '353939'}],
                                                    u'AssemblyAccession': 'GCA_001886755.1',
                                                    u'LastMajorReleaseAccession': 'GCA_001886755.1',
                                                    u'Synonym': {u'RefSeq': '',
                                                                 u'Genbank': 'GCA_001886755.1',
                                                                 u'Similarity': ''},
                                                    u'FtpPath_RefSeq': '',
                                                    u'AsmReleaseDate_RefSeq': '1/01/01 00:00',
                                                    u'AssemblyType': 'haploid',
                                                    u'AssemblyDescription': '',
                                                    u'AsmReleaseDate_GenBank': '2016/11/28 00:00',
                                                    u'RS_Projects': [],
                                                    u'RefSeq_category': 'na',
                                                    u'PartialGenomeRepresentation': 'false',
                                                    u'Coverage': '0',
                                                    u'AssemblyClass': 'haploid',
                                                    u'ExclFromRefSeq': [],
                                                    u'SeqReleaseDate': '2016/11/28 00:00',
                                                    u'AnomalousList': [],
                                                    u'AssemblyStatus': 'Complete Genome',
                                                    u'AsmReleaseDate': '2016/11/28 00:00',
                                                    u'ReleaseLevel': 'Major',
                                                    u'Taxid': '562',
                                                    u'LastUpdateDate': '2016/11/28 00:00',
                                                    u'RsUid': '',
                                                    u'FromType': '',
                                                    u'WGS': '',
                                                    u'GB_Projects': ['353939'],
                                                    u'BioSampleId': '6029894',
                                                    u'AssemblyName': 'ASM188675v1',
                                                    u'EnsemblName': '',
                                                    u'Organism': 'Escherichia coli (E. coli)',
                                                    u'Biosource': {u'Isolate': '',
                                                                   u'InfraspeciesList': [{u'Sub_type': 'strain',
                                                                                          u'Sub_value': 'MRSN346595'}],
                                                                   u'Sex': ''},
                                                    u'SpeciesName': 'Escherichia coli',
                                                    u'BioSampleAccn': 'SAMN06029894',
                                                    u'FtpPath_GenBank': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/886/755/GCA_001886755.1_ASM188675v1',
                                                    u'SpeciesTaxid': '562',
                                                    u'UCSCName': '',
                                                    u'Primary': '3745268',
                                                    u'Meta': ' <Stats> <Stat category="alt_loci_count" sequence_tag="all">0</Stat> <Stat category="chromosome_count" sequence_tag="all">1</Stat> <Stat category="contig_count" sequence_tag="all">6</Stat> <Stat category="contig_l50" sequence_tag="all">1</Stat> <Stat category="contig_n50" sequence_tag="all">4796421</Stat> <Stat category="non_chromosome_replicon_count" sequence_tag="all">5</Stat> <Stat category="replicon_count" sequence_tag="all">6</Stat> <Stat category="scaffold_count" sequence_tag="all">6</Stat> <Stat category="scaffold_count" sequence_tag="placed">6</Stat> <Stat category="scaffold_count" sequence_tag="unlocalized">0</Stat> <Stat category="scaffold_count" sequence_tag="unplaced">0</Stat> <Stat category="scaffold_l50" sequence_tag="all">1</Stat> <Stat category="scaffold_n50" sequence_tag="all">4796421</Stat> <Stat category="total_length" sequence_tag="all">5116627</Stat> <Stat category="ungapped_length" sequence_tag="all">5116627</Stat> </Stats> <FtpSites>   <FtpPath type="GenBank">ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/886/755/GCA_001886755.1_ASM188675v1</FtpPath> </FtpSites> <assembly-level>5</assembly-level> <assembly-status>Complete Genome</assembly-status> <representative-status>na</representative-status> <submitter-organization>Walter Reed Army Institute of Research</submitter-organization>    ',
                                                    u'NCBIReleaseDate': '2016/11/28 00:00',
                                                    u'SubmitterOrganization': 'Walter Reed Army Institute of Research',
                                                    u'RS_BioProjects': [],
                                                    u'SortOrder': '5C10018867559899'},
                                                    attributes={u'uid': u'894571'}
                                                   )
                                       ]
                  },
                 attributes={u'status': u'OK'}
                 )
 }
"""
