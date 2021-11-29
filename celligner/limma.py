###########################################################
#
# limmapy
#
##################################################################

from __future__ import print_function
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import rpy2
from rpy2.robjects.packages import importr
limma = importr('limma')
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro
import sys

to_dataframe = robjects.r('function(x) data.frame(x)')


class limmapy:
    '''
    limma object through rpy2
    input:
    count_matrix: should be a pandas dataframe with each column as count, and a id column for gene id
        example:
        id    sampleA    sampleB
        geneA    5         1
        geneB    4         5
        geneC    1         2
    design_matrix: an design matrix in the form of pandas dataframe, see limma manual, samplenames as rownames
                            treatment1, treatment2, ...
    sampleA        A          A          B
    sampleA        A          A          B
    sampleB        B          B          A
    sampleB        B          A          B
    '''

    def __init__(self):
        print("you need to have R installed with the limma library installed")
        print(rpy2.__version__)
        self.limma_result = None


    def lmFit(self, count_matrix, design_matrix, **kwargs):  # OPTIONAL
        """
        args:
            geoMeans: cond*gene matrix
        """
        with localconverter(ro.default_converter + pandas2ri.converter):
            count_matrix = pandas2ri.py2rpy(count_matrix.astype(int))
            design_matrix = pandas2ri.py2rpy(design_matrix.astype(int))
        self.fit = limma.lmFit(count_matrix, design_matrix, **kwargs)
        return self

    def eBayes(self, **kwargs):
        self.fit = limma.eBayes(self.fit, **kwargs)
        return self

    def topTable(self, **kwargs):
        val = limma.topTable(self.fit, **kwargs)
        if type(val) == robjects.vectors.DataFrame:
            with robjects.conversion.localconverter(
                    robjects.default_converter + pandas2ri.converter):
                val = ro.conversion.rpy2py(val)
        return val
