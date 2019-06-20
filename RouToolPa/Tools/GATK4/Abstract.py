#!/usr/bin/env python

__author__ = 'mahajrod'

from RouToolPa.Tools.Abstract import Tool


class GATKTool(Tool):

    def __init__(self, java_path="", max_threads=4, jar_path="", max_memory="1g"):

        jar = "GenomeAnalysisTK.jar"

        Tool.__init__(self, "java", path=java_path, max_threads=max_threads, jar_path=jar_path, jar=jar, max_memory=max_memory)
