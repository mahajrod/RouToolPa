__author__ = 'mahajrod'

from RouToolPa.Tools.Annotation.Emapper import Emapper
from RouToolPa.Tools.Annotation.VEP import VEP
from RouToolPa.Tools.Annotation.BUSCO import BUSCO
from RouToolPa.Tools.Annotation.SNPeff import SNPeff
from RouToolPa.Tools.Annotation.Barrnap import Barrnap
from RouToolPa.Tools.Annotation.AUGUSTUS import AUGUSTUS
from RouToolPa.Tools.Annotation.Exonerate import Exonerate
from RouToolPa.Tools.Annotation.TransDecoder import TransDecoder


SNPeff_path = ""
SNPeff = SNPeff(jar_path=SNPeff_path)

VEP = VEP()
BUSCO = BUSCO()
Barrnap = Barrnap()
Emapper = Emapper()
AUGUSTUS = AUGUSTUS()
Exonerate = Exonerate()
TransDecoder = TransDecoder()
