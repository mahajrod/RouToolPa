__author__ = 'mahajrod'

from RouToolPa.Tools.Alignment.EMA import EMA
from RouToolPa.Tools.Alignment.GEM import GEM
from RouToolPa.Tools.Alignment.BWA import BWA
from RouToolPa.Tools.Alignment.BLAT import BLAT
from RouToolPa.Tools.Alignment.TMAP import TMAP
from RouToolPa.Tools.Alignment.STAR import STAR
from RouToolPa.Tools.Alignment.GMAP import GMAP
from RouToolPa.Tools.Alignment.Tophat import Tophat
from RouToolPa.Tools.Alignment.BamUtil import BamUtil
from RouToolPa.Tools.Alignment.Bowtie2 import Bowtie2
from RouToolPa.Tools.Alignment.Mosdepth import Mosdepth
from RouToolPa.Tools.Alignment.Novoalign import Novoalign
from RouToolPa.Tools.Alignment.LongRanger import LongRanger

EMA = EMA()
GEM = GEM()
BWA = BWA()
TMAP = TMAP()
BLAT = BLAT()
GMAP = GMAP()
STAR = STAR()
BamUtil = BamUtil()
Tophat = Tophat()
Bowtie2 = Bowtie2()
Mosdepth = Mosdepth()
Novoalign = Novoalign()
LongRanger = LongRanger()
