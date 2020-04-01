from __future__ import division
from __future__ import unicode_literals
from IPython.core.display import display, HTML
from json import dumps as json_encode
from base64 import b64encode
from IPython.display import display, Javascript, clear_output
from IPython.core.debugger import set_trace
from subprocess import call
import os


# the base64 library produces clean output without newlines

call(['pip', 'install', 'allensdk'])
clear_output()
call(['pip', 'install', 'joblib'])
clear_output()
call(['pip', 'install', '--upgrade','hbp-service-client==1.0.0']) #!pip install --upgrade "hbp-service-client==1.0.0"
clear_output()
#!conda install -c r r-base -y
call(['pip', 'install', '--upgrade','sklearn']) #!pip install --upgrade sklearn
clear_output()
call(['git', 'clone', 'https://github.com/AllenInstitute/mouse_connectivity_models.git'])
clear_output()
os.chdir('mouse_connectivity_models/')
#call(['cd', 'mouse_connectivity_models/']) #!pip install --upgrade sklearn
clear_output()
call(['pip', 'install', '.']) #!pip install --upgrade sklearn
os.chdir('..')
clear_output()
call(['pip', 'install', 'pandas==0.19.2'])
clear_output()
call(['pip', 'install', 'nilearn'])
clear_output()
call(['pip', 'install', 'wget'])
clear_output()
import sys
if 'rpy2' not in sys.modules:
    try:
        import rpy2
	if rpy2.__version__!='2.8.6':
	   call(['pip', 'install', 'rpy2=2.8.6'])
	   import rpy2
    except:
        call(['pip', 'install', 'rpy2=2.8.6'])
        import rpy2
else:
	import rpy2
	if rpy2.__version__!='2.8.6':
	   call(['pip', 'install', 'rpy2=2.8.6'])
	   import rpy2

if 'numpy' not in sys.modules: 
    call(['pip', 'install', 'numpy=1.11.2'])       
    import numpy as np
else:
    import numpy as np
    if np.__version__ != '1.11.2':
	call(['pip', 'install', 'numpy=1.11.2'])       
	import numpy as np

import copy
import matplotlib
import scipy as sci
import csv
import pickle as pk
import h5py
import math as mt
import matplotlib.pyplot as plt
import matplotlib.axis as axis
import matplotlib
import matplotlib.legend as lgd
import sklearn.cluster as cluster
import seaborn as sns
import glob
import six
import time as time
import json
import pandas as pd
import mcmodels.core as core
import mcmodels.models as models
import sklearn.decomposition as decomp
import sklearn.feature_selection as fs
import sklearn.linear_model as lm
import sklearn.metrics as metrics
import zipfile
import nrrd

import rpy2.robjects as ro
from rpy2.robjects import numpy2ri,pandas2ri
import nibabel as nib
import base64
import skimage as ski
import imp
import io, sys, types
import xlrd as xlrd
import xlrd
import string
import urllib
import nilearn
import wget, gzip


from sklearn.decomposition import DictionaryLearning
from IPython import get_ipython
from nbformat import read
from IPython.core.interactiveshell import InteractiveShell
from sklearn import multioutput
from pandas.tools.plotting import table
from itertools import groupby
from sklearn.preprocessing import Imputer, normalize, StandardScaler, scale
from scipy.stats import boxcox, mstats, normaltest, pearsonr, skew, iqr
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, GradientBoostingClassifier
from sklearn.model_selection import cross_val_predict,\
                                    GridSearchCV, KFold,\
                                    permutation_test_score, \
                                    StratifiedKFold, cross_val_score, RandomizedSearchCV
from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.linear_model import RandomizedLogisticRegression,\
                                 LogisticRegression,SGDClassifier, Ridge, Lars, ElasticNet
from sklearn.neural_network import MLPRegressor
from sets import Set
from sklearn import utils
from collections import OrderedDict
from joblib import Parallel, delayed
from scipy.ndimage.filters import laplace
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from allensdk.api.queries.grid_data_api import GridDataApi
from allensdk.api.queries.ontologies_api import OntologiesApi
from mcmodels.models.voxel import RegionalizedModel
from subprocess import call
from itertools import chain
from os import chdir
from os.path import basename,dirname,realpath
from scipy.spatial.distance import euclidean
from scipy.stats import skew, kurtosis, zscore
from scipy.stats.mstats import normaltest
from allensdk.api.queries.grid_data_api import *
from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi
from allensdk.api.queries.reference_space_api import ReferenceSpaceApi
from scipy.io import loadmat
from urllib import urlretrieve
from sklearn.preprocessing import Imputer, normalize
from scipy.special import expit
from mcmodels.core import VoxelModelCache
from nilearn.image import load_img,resample_to_img
from sklearn.metrics import mean_squared_error,r2_score

os.environ['KMP_DUPLICATE_LIB_OK']='True'

import cortical_map_10 as cm
from NotebookSearch import *
