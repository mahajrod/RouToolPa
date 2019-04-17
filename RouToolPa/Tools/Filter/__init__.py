__author__ = 'mahajrod'

from RouToolPa.Tools.Filter.FaCut import FaCut
from RouToolPa.Tools.Filter.FastQC import FastQC
from RouToolPa.Tools.Filter.Adapters import adapters_PE
from RouToolPa.Tools.Filter.Cutadapt import Cutadapt
from RouToolPa.Tools.Filter.Clumpify import Clumpify
from RouToolPa.Tools.Filter.TrimGalore import TrimGalore
from RouToolPa.Tools.Filter.Trimmomatic import Trimmomatic
from RouToolPa.Tools.Filter.Cookiecutter import Cookiecutter
from RouToolPa.Tools.Filter.Cookiecutter import CookiecutterOld



FaCut = FaCut()
FastQC = FastQC()
Cutadapt = Cutadapt()
Clumpify = Clumpify()
TrimGalore = TrimGalore()
Trimmomatic = Trimmomatic()
Cookiecutter = Cookiecutter()
CookiecutterOld = CookiecutterOld()

