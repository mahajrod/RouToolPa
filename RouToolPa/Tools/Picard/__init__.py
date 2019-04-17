#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
from RouToolPa.Tools.Picard.MarkDuplicates import MarkDuplicates
from RouToolPa.Tools.Picard.SortVcf import SortVcf
from RouToolPa.Tools.Picard.AddOrReplaceReadGroups import AddOrReplaceReadGroups
from RouToolPa.Tools.Picard.CreateSequenceDictionary import CreateSequenceDictionary


SortVcf = SortVcf()
MarkDuplicates = MarkDuplicates()
AddOrReplaceReadGroups = AddOrReplaceReadGroups()
CreateSequenceDictionary = CreateSequenceDictionary()
