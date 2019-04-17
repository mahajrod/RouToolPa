__author__ = 'mahajrod'

from RouToolPa.Tools.Alignment.BLAT import BLAT
from RouToolPa.Tools.Alignment.GEM import GEM
from RouToolPa.Tools.Alignment.BWA import BWA
from RouToolPa.Tools.Alignment.TMAP import TMAP
from RouToolPa.Tools.Alignment.STAR import STAR
from RouToolPa.Tools.Alignment.Tophat import Tophat
from RouToolPa.Tools.Alignment.Bowtie2 import Bowtie2
from RouToolPa.Tools.Alignment.Novoalign import Novoalign
from RouToolPa.Tools.Alignment.LongRanger import LongRanger


max_threads = 4
bowtie2_path = ""
bwa_path = ""
novoalign_path = ""
tmap_path = ""
prank_path = ""

GEM = GEM()
BWA = BWA(path=bwa_path, max_threads=max_threads)
TMAP = TMAP(path=novoalign_path, max_threads=max_threads)
BLAT = BLAT()
STAR = STAR()
Tophat = Tophat()
Bowtie2 = Bowtie2(path=bowtie2_path, max_threads=max_threads)
Novoalign = Novoalign(path=novoalign_path, max_threads=max_threads)
LongRanger = LongRanger()
