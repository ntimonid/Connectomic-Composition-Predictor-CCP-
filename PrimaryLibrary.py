import sys
import os
import numpy as np
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
import cortical_map_10 as cm
import nrrd
import rpy2.robjects as ro
import rpy2
import nibabel as nib
import base64
import nrrd
import skimage as ski

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
from scipy.ndimage.filters import laplace
from allensdk.api.queries.ontologies_api import OntologiesApi
from IPython.core.display import display, HTML
from __future__ import unicode_literals
from IPython.display import display, Javascript, clear_output
from json import dumps as json_encode
from base64 import b64encode
from mcmodels.models.voxel import RegionalizedModel
from IPython.core.debugger import set_trace

os.environ['KMP_DUPLICATE_LIB_OK']='True'

class MesoconnectomePredictor:
    
    __slots__ = ['GeneExp','ConStr']
    def __init__(self, GeneExp = [], ConStr = []):
        self.GeneExp = GeneExp
        self.ConStr = ConStr
        
        
    def GetLayerResolvedArray(cre_file = 'Supplementary Table 1.csv', creFilter = True):
        
        if creFilter == True:
            # Default filtering of cre-lines based on the 15 tracer lines used by Harris et al 2018
            creFilter = ['Syt6-Cre_KI148', 'Ntsr1-Cre_GN220', 'Sim1-Cre_KJ18',
                        'Efr3a-Cre_NO108', 'Chrna2-Cre_OE25', 'A93-Tg1-Cre',  
                        'Tlx3-Cre_PL56', 'Rbp4-Cre_KL100', 'Rorb-IRES2-Cre',
                        'Scnn1a-Tg3-Cre', 'Nr5a1-Cre', 'Sepw1-Cre_NP39',
                        'C57BL/6J', 'Emx1-IRES-Cre', 'Cux2-IRES-Cre'] 


        creDict = {}
        with open(cre_file) as fp:
            buff = csv.reader(fp)
            for idx,row in enumerate(buff):
                 if len(row) > 1: # concatenate the two rows together - error caused by csv transition
                    row[0] =  row[0] + row[1]
                 remains = [x for x in filter(None, row[0].split(';'))]
                 if remains[0].isdigit():
                    if remains[1] in creFilter or creFilter == False:
                        if 'A93' in remains[1]:         # name incoherence between the default creFilter and the cre_file for A93
                            remains[1] = 'A930038C07Rik-Tg1-Cre'
                        if '-' in remains[7]:           #  covers layers 2-6 and therefore is layer inspecific
                            profile = 'layer inspecific'
                        else:    
                            profile = remains[7] + ' ' + remains[8]
                        if profile not in creDict:
                            creDict[profile] = []
                            creDict[profile].append(remains[1])
                        else:                
                            creDict[profile].append(remains[1])

        cache = core.VoxelModelCache()
        layer_resolved_array = []

        for cur_profile in creDict.keys():
            print 'Current analyzed cre lines are: {} -> {}'.format(cur_profile, creDict[cur_profile])

            # gather injection and projection data from all tracing experiments
            # belonging to the selected cre-line
            tracing_data = cache.get_experiment_data(cre = creDict[cur_profile],
                                                     injection_hemisphere_id  = 3,
                                                     projection_hemisphere_id = 3,
                                                     flip_experiments = True,
                                                     injection_structure_ids  = [315])

            source_voxels  = tracing_data.injection_mask.coordinates
            injections_64  = np.asarray(tracing_data.injections, dtype = np.float64)
            projections_64 = np.asarray(tracing_data.projections, dtype = np.float64)
            source_key = tracing_data.injection_mask.get_key()
            target_key = tracing_data.projection_mask.get_key() 
            # create a voxel-scale model of the mouse connectome given the injection
            # and projection data
            voxelModel    = models.VoxelModel(source_voxels)
            voxelModel.fit((tracing_data.centroids, injections_64), projections_64)

            # generate a voxel-scale connectivity array based on the fit model
            voxel_array = models.voxel.VoxelConnectivityArray.from_fitted_voxel_model(voxelModel)

            # regionalize model
            regional_model = models.voxel.RegionalizedModel.from_voxel_array(
                 voxel_array, source_key, target_key)

            layer_resolved_array.append(regional_model.normalized_connection_density.T)
            
        # creation of a 3D array serving our laminarly resolved connectivity array    
        layer_resolved_array = np.asarray(layer_resolved_array, dtype = np.float32)

        return layer_resolved_array    
             
        
    def GOenrichment(self, GeneList):
        # Description: given a set of selected genes, performs enrichment
        #              analysis comparing to a bioconductor database
        ro.numpy2ri.activate()
        MGP = ro.r('''
        GOenrichment <-function(GeneList){
            source("http://bioconductor.org/biocLite.R")
            biocLite("ALL")
            biocLite("GOstats")
            biocLite("Category")
            biocLite("genefilter")
            biocLite("org.Mm.eg.db")
            biocLite("gage")
            biocLite(pkgs=c("Biobase", "IRanges", "AnnotationDbi"),
                suppressUpdates=FALSE,
                suppressAutoUpdate=FALSE,
                siteRepos=character(),
                ask=TRUE)
            #install.packages("RSQLite")
            #install.packages("devtools")
            #require(devtools)
            library("RSQLite")
            library("gage")
            library("genefilter")
            library("org.Mm.eg.db")
            library(DBI)
            library("GO.db")
            library("GOstats")
            library("Category")
            library("annotate")
            data(ALL, package="ALL")

            hgCutoff <- 0.001
            params <- new("GOHyperGParams",
                   geneIds = GeneList,
                   universeGeneIds = NULL,
                   annotation = "org.Mm.eg.db",
                   ontology = "BP",
                   pvalueCutoff = hgCutoff,
                   conditional = FALSE,
                   testDirection = "over")
            BP  <- hyperGTest(params)
            sumBP <- data.frame(summary(BP))
            sumBP <- subset(sumBP, select=c("Term"))

            params <- new("GOHyperGParams",
                   geneIds = GeneList,
                   universeGeneIds = NULL,
                   annotation = "org.Mm.eg.db",
                   ontology = "MF",
                   pvalueCutoff = hgCutoff,
                   conditional = FALSE,
                   testDirection = "over")
            MF <- hyperGTest(params)
            sumMF <- data.frame(summary(MF))
            sumMF <- subset(sumMF, select=c("Term"))

            params <- new("GOHyperGParams",
                   geneIds = GeneList,
                   universeGeneIds = NULL,
                   annotation = "org.Mm.eg.db",
                   ontology = "CC",
                   pvalueCutoff = hgCutoff,
                   conditional = FALSE,
                   testDirection = "over")
            CC <- hyperGTest(params)
            sumCC <- data.frame(summary(CC))
            sumCC <- subset(sumCC, select=c("Term"))

            foo <- vector(mode="list", length = 3)
            foo[[1]] <- sumBP
            foo[[2]] <- sumMF
            foo[[3]] <- sumCC

            return(foo)
        }''')
        r_getname = ro.globalenv['GOenrichment']
        verdict = r_getname(GeneList)
        BP = verdict[0]
        MF = verdict[1]
        CC = verdict[2]

        return BP, MF, CC    
        
        
        
    def PlotTheResults(self,Results,MetaInfo):

            k = len(MetaInfo['xtick'])
            fsz = 9
            fig = plt.figure(figsize=(8,6))

            #plt.tight_layout()
            plt.rcParams['figure.figsize']
            plt.boxplot(Results, 0, 'gD', widths = 0.6, whis = [5,95])

            ax = plt.gca()
            ax.xaxis.set_tick_params(labelsize = 14)
            ax.yaxis.set_tick_params(labelsize = 14)
            plt.xticks([i+1 for i in range(k)], MetaInfo['xtick'],rotation = 0)
            plt.text(0.55, MetaInfo['max1'][1], MetaInfo['max1'][0], fontsize = fsz)
            plt.text(0.55, MetaInfo['min1'][1], MetaInfo['min1'][0], fontsize = fsz)
            plt.text(0.55, MetaInfo['med1'][1], MetaInfo['med1'][0], fontsize = fsz)
            plt.text(1.45, MetaInfo['max2'][1], MetaInfo['max2'][0], fontsize = fsz)
            plt.text(1.45, MetaInfo['min2'][1], MetaInfo['min2'][0], fontsize = fsz)
            plt.text(1.45, MetaInfo['med2'][1], MetaInfo['med2'][0], fontsize = fsz)
            plt.text(2.45, MetaInfo['max3'][1], MetaInfo['max3'][0], fontsize = fsz)
            plt.text(2.45, MetaInfo['min3'][1], MetaInfo['min3'][0], fontsize = fsz)
            plt.text(2.45, MetaInfo['med3'][1], MetaInfo['med3'][0], fontsize = fsz)
            plt.title(MetaInfo['title'], fontsize = 14)
            # plt.xlabel(MetaInfo['xlabel'])
            plt.ylabel(MetaInfo['ylabel'], fontsize = 14)
            plt.yticks(np.arange(MetaInfo['lb'] - np.std(Results), MetaInfo['ub'] + np.std(Results), np.std(Results)))
            #plt.yticks(np.arange(0.4,1.0,0.05))

            plt.savefig(MetaInfo['save_file'])
            plt.show(block = False)
            plt.pause(1)
            plt.close() 
            
             
    def Evaluation(self,ClfResults_lr,ClfResults_rf, ClfResults_bl, params):
 
          MetaInfo = {}
          templ = 'values for predicted tracer - projection patterns\nof driver line:'
          MetaInfo['xlabel'] = 'Classifier Models'
          exclusions = ['Model Storage', 'Gene Scoring', 'Mean Gene Scoring']
          for measure in ClfResults_lr.keys():
             if measure not in exclusions:
               tmp1 = [val[0] for val in ClfResults_lr[measure]] 
               tmp2 = [val[0] for val in ClfResults_rf[measure]] 
               tmp3 = [val[0] for val in ClfResults_bl[measure]] 
               sort1 = np.argsort(tmp1)[::-1]; sort3 = np.sort(tmp1)[::-1]
               sort2 = np.argsort(tmp2)[::-1]; sort4 = np.sort(tmp2)[::-1]
               sort5 = np.argsort(tmp3)[::-1]; sort6 = np.sort(tmp3)[::-1]
               set_length = len(sort1)   
            
               MetaInfo['min1']  = [params['structure-abbrev'][sort1[set_length-1]], sort3[set_length-1]]
               MetaInfo['max1']  = [params['structure-abbrev'][sort1[0]], sort3[0]]
               MetaInfo['med1']  = [params['structure-abbrev'][sort1[set_length/2]], sort3[set_length/2]]
    
               MetaInfo['min2']  = [params['structure-abbrev'][sort2[set_length-1]], sort4[set_length-1]]
               MetaInfo['max2']  = [params['structure-abbrev'][sort2[0]], sort4[0]]
               MetaInfo['med2']  = [params['structure-abbrev'][sort2[set_length/2]], sort4[set_length/2]]
               MetaInfo['min3']  = [params['structure-abbrev'][sort5[set_length-1]], sort6[set_length-1]]
               MetaInfo['max3']  = [params['structure-abbrev'][sort5[0]], sort6[0]]
               MetaInfo['med3']  = [params['structure-abbrev'][sort5[set_length/2]], sort6[set_length/2]] 
               
               MetaInfo['save_file'] = params['cre-line'] + '_' + measure +'.png'
               MetaInfo['title']     = params['cre-line'] #measure + ' '+ templ + ' ' + params['cre-line']
               MetaInfo['ylabel']    = measure
               MetaInfo['xtick']     = ['Ridge Regression','Random Forest', 'Control']
               if measure == 'AURoc': 
                  MetaInfo['lb'] = 0.4
                  MetaInfo['ub'] = 1.0
               elif measure == 'RMSE' or measure == 'r2': 
                  MetaInfo['lb'] = min([sort3[set_length-1],sort4[set_length-1],sort6[set_length-1]])
                  MetaInfo['ub'] = max([sort3[0],sort4[0],sort6[0]])
               else: 
                  MetaInfo['lb'] = 0.1 
                  MetaInfo['ub'] = 1.0
                
               JoinResults = np.asarray([(a[0],b[0],c[0]) for a,b,c in zip(ClfResults_lr[measure],ClfResults_rf[measure],ClfResults_bl[measure])])
               MesoPred.PlotTheResults(JoinResults,MetaInfo)  
          
          '''gene_stars =  50
          BP, MF, CC = MesoPred.GOenrichment(ClfResults_lr['Mean Gene Scoring'][0][0:gene_stars])
          GODict = OrderedDict()
          GODict['BIOLOGICAL PROCESS'] = [element for element in BP]
          GODict['MOLECULAR FUNCTION'] = [element for element in MF]
          GODict['CELLULAR COMPONENT'] = [element for element in CC]
          GO_DF = pd.DataFrame(dict([ (k, pd.Series(v)) for k,v in GODict.items() ])) 
          table(ax, GO_DF)'''
          
        
        
        
    def StabilityInvestigator(model_list, unanimity_thr):
    
        frequency = {}; stability_set = []; params_to_fit = []
        param_dict = [model.best_params_ for model in model_list]
        param_list = [fold.values() for fold in param_dict]

        for param in set(chain.from_iterable(param_list)):
            frequency[param] = len([idx for idx,model_params in enumerate(param_list) if param in model_params])
            if frequency[param] == len(param_dict): stability_set.append(param)

        if len(stability_set)/(1.0*len(np.unique(param_list))) > unanimity_thr:
            print 'model is stable, train the whole dataset with the parameters'
            params_to_fit = param_dict[0]
            verdict = 'Full Fit'
        elif len(stability_set)/(1.0*len(np.unique(param_list))) > unanimity_thr - 0.2:
            verdict = 'Emsemble Fit'
            print 'model is slightly unstable, ensemble training over the cross-validation models'
        else:
            verdict = 'No Fit'
            print 'model is highly unstable, possibly by overfitting'

        return verdict, params_to_fit
    
    
    def PreProcessing(self,GeneExp, ConStr, params):
        
        # ************ Step 1: Cleaning the dataset from rows containing complete nans *******************#
        print 'Data cleaning phase ...\n'
        print 'Size before data cleaning ' + str(np.shape(ConStr))
        
        ConStr,nanCons = RemoveNanStructs(0,0.1).fit(ConStr)
        print 'Size after intermediate cleaning ' + str(np.shape(ConStr))
        params['restCons'] = [idx for idx in range(len(GeneExp)) if idx not in nanCons]
        GeneExp = np.delete(GeneExp,nanCons,0)
       
        GeneExp,nanGenes = RemoveNanStructs(0,0.03).fit(GeneExp)
        params['nanGenes'] = nanGenes; params['nanCons'] = nanCons
        params['str_cpy'] = np.delete(params['structure-abbrev'],params['nanCons'],0); params['str_cpy'] = np.delete(params['str_cpy'],params['nanGenes'],0)
        params['restGenes'] = [idx for idx in range(len(ConStr)) if idx not in nanGenes]
        ConStr = np.delete(ConStr,nanGenes,0)
        print 'Size after data cleaning ' + str(np.shape(ConStr))
        m = len(ConStr[0])
        tmp = [params['leaf_keys'][val] for val in params['restCons']]
        params['remaining_indices']   = [tmp[val] for val in params['restGenes']]
 
        # ************** Impute the value of the remaining nans with column averaging **************#
        print 'Data imputation phase ...\n'
        for task in range(len(ConStr[0])):
            nans = [val for val in ConStr[:,task] if mt.isnan(val) == True]
            if len(nans) > 0:
                ConStr[:,task] = BimodalImputation().fit(ConStr[:,task])

        GeneExp = Imputer(missing_values = 'NaN', strategy = 'median', axis = 0, verbose=0, copy=True).fit_transform(GeneExp)

        print 'Data normalization phase ....\n'
        
        GeneExp_hat = np.power(GeneExp, 1./3)
        GeneExp_hat_ctr = StandardScaler().fit_transform(GeneExp_hat)
        ConStr_sqrt = np.power(ConStr, 1./3)
        sc_sqrt = StandardScaler()
        sc_sqrt.fit(ConStr_sqrt)
        ConStr_sqrt_ctr = sc_sqrt.transform(ConStr_sqrt)

        # Remove outliers
        GeneExp_hat_ctr,rem_genes = RemoveOutliers().fit(GeneExp_hat_ctr)
        params['Gene Acronyms'] = np.delete(params['Gene Acronyms Original'],rem_genes,0);
        params['rem_genes'] = rem_genes
        print 'Removing outliers ... Reduced dim: '+ str(np.shape(GeneExp)) + ' ' + str(np.shape(ConStr))


        # Binarize the connectivity matrix for classification
        binThr  = 0.33333
        y       = np.zeros((np.shape(ConStr)))
        for j in range(ConStr.shape[1]):
            y[:,j] = BinarizeTheVector(ConStr[:,j],binThr)

        return GeneExp_hat_ctr, ConStr_sqrt_ctr, y, sc_sqrt
        
        
        
    def LaminarRegistration(acronyms):
 
        laminar_profiles = []
        for acro in acronyms:
            numbers = [idx for idx,val in enumerate(acro) if val.isdigit()]
            if numbers != []: 
                split_letter =  acro[numbers[0]-1]
                split_choices =  acro.split(split_letter)
                potential_choice = split_choices[len(split_choices)-1]
                blobs = [val for val in potential_choice if val.isdigit() == False]
                if len(blobs) < 2 and potential_choice != '' and 'r' not in potential_choice: 
                    laminar_profiles.append('layer ' + potential_choice)
                elif '6a' in acro:
                    laminar_profiles.append('layer 6a')
                else:
                    laminar_profiles.append('layer inspecific')
            else:
                laminar_profiles.append('layer inspecific')

        return laminar_profiles   
    
    def Convert2ROC(self,t_actual,t_pred):
        
        def binarize(val,Q):
            if val < Q:
                return 0
            else:
                return 1
        
        cutoff_collector = {}
        cutoffs = np.floor(np.linspace(0,100,100))
        auc_vec = np.zeros((np.shape(t_actual)[1], len(cutoffs)))
        for tracer in range(np.shape(t_actual)[1]):
            for cutoff in cutoffs:
                Q = np.percentile(list(t_actual[:,tracer]), 100-cutoff, axis = 0)
                y_actual = np.vectorize(binarize)(t_actual[:,tracer],Q)
                y_pred   = 1/(1 + np.exp(- t_pred[:,tracer]))
                ones = [val for val in y_actual if val == 1]
                if cutoff >= 20 and cutoff <= 80 and len(np.unique(y_actual)) == 2:
                    auc = metrics.roc_auc_score(y_actual, y_pred, average = 'micro')
                    auc_vec[tracer,int(cutoff)] = auc
             
            cutoff_collector[tracer] = (cutoffs[np.argmax(auc_vec[tracer,:])], np.max(auc_vec[tracer,:]))    

        return cutoff_collector          
        
            
    def UnravelResults(self,predictions, params, scaler):
         
        num      = 1
        if len(predictions) == 1:
            predictions = predictions[0]
        if params['method'] == 'classification':
            y        = np.transpose([task[2] for task in predictions])
            y_preds  = np.transpose([task[1] for task in predictions])
            y_scores = np.asarray([task[0] for task in predictions])
            y_scores = y_scores.transpose(1,0,2)
        elif params['method'] == 'regression':
            if len(predictions) > 5:
                y        = np.transpose([task[2] for task in predictions])
                y_preds  = np.transpose([task[1] for task in predictions])
                y_scores = np.asarray([task[0] for task in predictions])
                if len(predictions[0]) > 4:
                    coeffs = [task[4] for task in predictions]
                mdls     = [task[3] for task in predictions]    
            else:
                y        = predictions[2]
                y_preds  = predictions[1]
                y_scores = predictions[0]
                mdls     = predictions[3]
                coeffs   = predictions[4]
                y_preds  = scaler.inverse_transform(y_preds)
                y        = scaler.inverse_transform(y)
                    
        if params['distance prior'] == True:
           y_preds  = np.asarray([MesoPred.DiscrimFunc(y_score, target, params) for target, y_score in enumerate(y_scores)]).transpose()

        ClfResults = {}
        ClfResults['Model Storage'] = mdls
        if params['method'] == 'classification':
            ClfResults['Precision'] = np.zeros((np.shape(y)[num],1))
            ClfResults['F-score']   = np.zeros((np.shape(y)[num],1))
            ClfResults['Recall']    = np.zeros((np.shape(y)[num],1))
            ClfResults['AURoc']     = np.zeros((np.shape(y)[num],1))
        elif params['method'] == 'regression':
            ClfResults['RMSE']      = np.zeros((np.shape(y)[num],1))
            ClfResults['r2']        = np.zeros((np.shape(y)[num],1))
      
        if len(predictions) > 4 and np.shape(coeffs)!= ():
            mean_coeffs = np.mean(coeffs,0)
            tmp_mean = np.argsort(mean_coeffs)[::-1]
            ClfResults['Mean Gene Scoring'] = [params['Gene Acronyms'][tmp_mean], mean_coeffs[tmp_mean], tmp_mean]
            ClfResults['Gene Scoring'] = []
            for coeff in coeffs:
                tmp = np.argsort(coeff)[::-1]
                ClfResults['Gene Scoring'].append([params['Gene Acronyms'][tmp],coeff[tmp], tmp])    
                
        '''elif params['method'] == 'regression' and len(predictions) > 4 and np.shape(coeffs)!= ():
            tmp = np.argsort(coeffs)[::-1]
            ClfResults['Gene Scoring'] = [params['Gene Acronyms'][tmp], coeffs[tmp], tmp]'''
       
        for num in range(len(y[0])):
            injection_task = 'Task ' + str(num)
            if params['method'] == 'classification':
                ClfResults['F-score'][num]   = metrics.f1_score(y[:,num], y_preds[:,num], \
                                                                average = 'binary')
                ClfResults['Precision'][num] = metrics.precision_score(y[:,num], \
                                                               y_preds[:,num],\
                                                               average = 'binary')
                ClfResults['Recall'][num]    = metrics.recall_score(y[:,num], y_preds[:,num], \
                                                             average = 'binary')
                ClfResults['AURoc'][num]     = metrics.roc_auc_score(y[:,num],y_scores[:,num][:,1],\
                                                             average = 'micro')
            elif params['method'] == 'regression':
                ClfResults['RMSE'][num]      = metrics.mean_squared_error(y[:,num], y_preds[:,num])
                ClfResults['r2'][num]        = metrics.r2_score(y[:,num], y_preds[:,num])

        return ClfResults    
    
    
    
    
    
    
    
    
    