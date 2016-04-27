__author__ = 'KATRINA'

'''

Run ZIFA algorithm on input file, print the result; this is intended to test the method.
Implementation of the result was done in ipython notebook.

ZIFA_test1.py [input_file]

'''

import pandas as pd
import numpy as np
import sys

from ZIFA import ZIFA
#from ZIFA import block_ZIFA

print('script started')

input_file = sys.argv[1]
df = pd.DataFrame.from_csv((input_file),sep="\t")
f = lambda x: np.log(1+x)
df1 = df.applymap(f)
print(df1)
print('completed read in DF')

#Z, model_params = block_ZIFA.fitModel(df1.as_matrix(), 2)
Z, model_params = ZIFA.fitModel(df1.as_matrix(), 2)

print('ZIFA finished')
print(Z)

#print(Z[:0])
#print(Z[:1])