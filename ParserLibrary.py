import numpy as np
import os
import pickle as pk
import h5py
import sys
import math as mt
import csv
import matplotlib.pyplot as plt
import numpy as np
import xlrd
import nrrd
import numpy as np
import json
import string
import csv
import scipy as sci
from IPython.core.debugger import set_trace
os.environ['KMP_DUPLICATE_LIB_OK']='True'
from os import chdir
from os.path import basename,dirname,realpath
from scipy.spatial.distance import euclidean
from scipy.stats import skew, kurtosis, zscore
from scipy.stats.mstats import normaltest
from allensdk.api.queries.grid_data_api import *
from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi
from allensdk.api.queries.reference_space_api import ReferenceSpaceApi
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from scipy.io import loadmat
from subprocess import call
from urllib import urlretrieve
from sklearn.preprocessing import Imputer, normalize
from collections import OrderedDict
from scipy.special import expit
import nrrd
from mcmodels.core import VoxelModelCache
from mcmodels.models.voxel import RegionalizedModel


API_PATH = "http://api.brain-map.org/api/v2/data"
GRAPH_ID = 1
PLANE_ID = 1 # coronal
MOUSE_PRODUCT_ID = 1 # aba

DATA_CON_SET_QUERY_URL = ("%s/SectionDataSet/query.json" +\
                          "?criteria=[failed$eqfalse]" +\
                          ",products[id$in5]" +\
                          ",[green_channel$eqrAAV]" +\
                          #",specimen(donor[transgenic_mouse_id$eqall donors])" +\
                          ",specimen(stereotaxic_injections(age[days$ge54],[days$le58]))" +\
                          ",plane_of_section[id$eq%d]" +\
                          "&include=specimen(stereotaxic_injections(age,stereotaxic_injection_materials,stereotaxic_injection_coordinates,primary_injection_structure)),specimen(donor(age))") \
                          % (API_PATH, PLANE_ID)

UNIONIZE_CON_FMT = "%s/ProjectionStructureUnionize/query.json" +\
               "?criteria=[section_data_set_id$eq%d],[is_injection$eqfalse]" +\
               "&include=hemisphere"

STRUCTURES_URL = ("%s/Structure/query.json?" +\
                      "criteria=[graph_id$eq%d]") \
                      % (API_PATH, GRAPH_ID)

    
DATA_EXP_SET_QUERY_URL = ("%s/SectionDataSet/query.json" +\
                          "?criteria=[failed$eq'false'][expression$eq'true']" +\
                          ",products[id$eq%d]" +\
                          ",plane_of_section[id$eq%d],genes" +\
                          "&include=genes") \
                          % (API_PATH, MOUSE_PRODUCT_ID, PLANE_ID)

UNIONIZE_EXP_FMT = "%s/StructureUnionize/query.json" +\
               "?criteria=[section_data_set_id$eq%d],structure[graph_id$eq1]" +\
               ("&include=section_data_set(products[id$in%d])" % (MOUSE_PRODUCT_ID)) +\
               "&only=id,structure_id,sum_pixels,expression_energy,section_data_set_id"
    
    
    
    
# Make a query to the API via a URL.
def QueryAPI(url):
    start_row = 0
    num_rows = 2000
    total_rows = -1
    rows = []
    done = False

    # The ontology has to be downloaded in pages, since the API will not return
    # more than 2000 rows at once.
    while not done:
        pagedUrl = url + '&start_row=%d&num_rows=%d' % (start_row,num_rows)

        print pagedUrl
        source = urllib.urlopen(pagedUrl).read()

        response = json.loads(source)
        rows += response['msg']

        if total_rows < 0:
            total_rows = int(response['total_rows'])

        start_row += len(response['msg'])

        if start_row >= total_rows:
            done = True

    print('Number of results: {}'.format(total_rows))
    return rows



# Download the mouse brain structures in a structure graph.
def DownloadStructures():
    structs = QueryAPI(STRUCTURES_URL)

    # Build a dict from structure id to structure and identify each node's
    # direct descendants.
    structHash = {}
    for s in structs:
        s['num_children'] = 0
        s['structure_id_path'] = [int(sid) for sid in s['structure_id_path'].split('/') if sid != '']
        structHash[s['id']] = s

    for sid,s in structHash.iteritems():
        if len(s['structure_id_path']) > 1:
            parentId = s['structure_id_path'][-2]
            structHash[parentId]['num_children'] += 1

    ## pull out the structure ids for structures in this structure graph that
    ## have no children (i.e. just the leaves)
    ## corrStructIds = [sid for sid,s in structHash.iteritems() if s['num_children'] == 0]
    # RB: no, leave all structures in and filter later
    corrStructIds = structHash.keys()

    return sorted(corrStructIds), structHash



def CreateConnectivityMatrix(dataSets,structureIds,structHash,unionizes):
    # Each injection experiment will have a connectivity vector.  This vector will be as long
    # as the number of requested structures.
    nstructs = len(structureIds)
    ndata = len(unionizes)
    print('ndata {} ndatasets {}'.format(ndata,len(dataSets)))

    sidHash = dict([(id,i) for (i,id) in enumerate(structureIds)])
    didHash = dict([(d['id'],i) for (i,d) in enumerate(dataSets)])

    connectivityL = numpy.empty([nstructs,ndata])
    connectivityL.fill(numpy.nan)
    connectivityR = numpy.empty([nstructs,ndata])
    connectivityR.fill(numpy.nan)
 
    connectivityDict = {'projection_density': 0, 'projection_intensity': 0, 'projection_energy': 0, 'projection_volume': 0, 'normalized_projection_volume': 0}
    for key in connectivityDict.keys():
        connectivityDict[key] = numpy.empty([nstructs,ndata])
        connectivityDict[key].fill(numpy.nan)

    # For each data set's set of unionizes, then for each individual structure,
    # fill in the structure's connectivity vector.
    for i,us in enumerate(unionizes):
        # for each unionize
        for j,u in enumerate(us):
            sid = u['structure_id']
            did = u['section_data_set_id']

            struct = structHash[sid]
            struct['volume'] = u['sum_pixels']

            if i ==0 and j == 0:
              print u

            if sidHash.has_key(sid) and didHash.has_key(did):
                if u['hemisphere_id'] is 1:
                    connectivityL[sidHash[sid]][didHash[did]]  = u['normalized_projection_volume']
                elif u['hemisphere_id'] is 2:
                    connectivityR[sidHash[sid]][didHash[did]] = u['normalized_projection_volume']
                    for key in connectivityDict.keys():
                        connectivityDict[key][sidHash[sid]][didHash[did]] = u[key]
                elif u['hemisphere_id'] is 3:
                  pass
                  # this is just the average value of L+R
            else:
                print "ERROR: structure {}/injection {} skipped.".format(sid,did)

    return connectivityL, connectivityR, connectivityDict

def CreateExpressionMatrix(dataSets,structureIds,structHash,unionizes):
    # Each structure will have an expression vector.  This vector will be as long
    # as the number of requested structures.
    nstructs = len(structureIds)
    ndata = len(unionizes)

    sidHash = dict([(id,i) for (i,id) in enumerate(structureIds)])
    didHash = dict([(d['id'],i) for (i,d) in enumerate(dataSets)])

    expression = numpy.empty([nstructs,ndata])
    expression.fill(numpy.nan)

    # For each data set's set of unionizes, then for each individual structure,
    # fill in the structure's expression vector.
    for i,us in enumerate(unionizes):
        # for each unionize
        for j,u in enumerate(us):
            sid = u['structure_id']
            did = u['section_data_set_id']

            struct = structHash[sid]
            struct['volume'] = u['sum_pixels']

            if sidHash.has_key(sid) and didHash.has_key(did):
                expression[sidHash[sid]][didHash[did]] = u['expression_energy']

    return expression



def allChildren(acr,acr2parent):
   # Description: given a tree hierarchy and an entity,
   # this function returns all the children
   # of the entity

   AC = []
   for a,p in acr2parent.items():
     if p == acr:
       AC.append(a)
       AC.extend(allChildren(a,acr2parent))
   return AC

def ReduceToLeafNodes(structure_acronyms,tree_file):
    # Description: this function checks the givenReduceToLeafNodes structures
    #              based on the tree hierarchy and returns the
    #              leaf nodes

    leaf_nodes = []
    with open(tree_file) as fp:
        acr2parent = json.load(fp)

    for idx,acro in enumerate(structure_acronyms):
        AC = allChildren(acro, acr2parent)
        if len(AC) == 0: # structure is a leaf node
           leaf_nodes.append((idx,acro))
    return leaf_nodes


def GetConUnionizes():

    infile1 = 'expression_files/inj_unionizes.nrrd'
    infile2 = 'expression_files/exp_density.nrrd'
    mcc = MouseConnectivityCache(resolution = 100)
    experiments = mcc.get_experiments(cre = True, dataframe=True)
    uni_con = []
    for idx,val in enumerate(experiments['id']):
        print idx
        tmp = mcc.get_experiment_structure_unionizes(experiment_id = val)
        uni_con.append(tmp)
    print uni_con.shape
    fp = h5py.File('py_files/unionized_connectivity.hdf5','w')
    fp.create_dataset('dataset1',data = uni_con)


def GetCreLines():
    infile = 'cre_inj_density.nrrd'
    infile2 = 'cre_pro_density.nrrd'
    mcc = MouseConnectivityCache(resolution = 100)
    mca = MouseConnectivityApi()
    cre_experiments2 = mca.experiment_source_search(injection_structures = 'root', transgenic_lines= True)
    cre_experiments = mcc.get_experiments(cre= True,dataframe=True)
    MetaPerCre = []
    InjPerCre  = []
    ProjPerCre = []

    creDict = {}
    with open('Supplementary Table 1.csv') as fp:
         buff = csv.reader(fp)
         for idx,row in enumerate(buff):
             if len(row) > 1: # concatenate the two rows together - error caused by csv transition
                row[0] =  row[0] + row[1]
             remains = [x for x in filter(None,row[0].split(';'))]
             if remains[0].isdigit():
             #if len(remains) > 8 and idx > 2 and remains[1].isdigit() == False:
                creDict[remains[1]] = []
                creDict[remains[1]].append(remains[7])
                creDict[remains[1]].append(remains[8])

    #pk.dump(cre_experiments,open('py_files/cre_experiments.pkl','wb'))
    #cre_experiments = pk.load(open('py_files/cre_experiments.pkl','rb'))

    # download the projection density volume for one of the experiments
    for idx,val in enumerate(cre_experiments['id']):
        mcc.get_projection_density(val, infile2)
        tmp = cre_experiments['transgenic-line'][val]
        if 'A93' in cre_experiments['transgenic-line'][val]:
            tmp = 'A93-Tg1-Cre'
        selCre = [key for key in creDict.keys() if tmp == key]
        if len(selCre) > 0:
            selCre = selCre[0]
            MetaPerCre.append({})
            rx = len(MetaPerCre)-1
            MetaPerCre[rx]['injection-coordinates'] = \
            cre_experiments['injection-coordinates'][val]
            MetaPerCre[rx]['structure-abbrev'] =\
             cre_experiments['structure-abbrev'][val]
            MetaPerCre[rx]['transgenic-line'] = tmp
            MetaPerCre[rx]['id'] = cre_experiments['id'][val]
            MetaPerCre[rx]['layer'] = creDict[selCre][0]
            MetaPerCre[rx]['Cell Type'] = creDict[selCre][1]
            # read it into memory
            pd_array, pd_info = nrrd.read(infile2)
            ProjPerCre.append(pd_array)

    f2 = h5py.File('ProjPerCre.hdf5','w')
    f2.create_dataset('dataset1',data = ProjPerCre)
    pk.dump(MetaPerCre,open('MetaPerCre.pkl','wb'))
    return MetaPerCre, ProjPerCre



def ReadConnectivityData():
    infile = 'expression_files/inj_density.nrrd'
    infile2 = 'expression_files/proj_density.nrrd'
    mca = MouseConnectivityApi()
    mcc = MouseConnectivityCache(resolution = 100)
    all_experiments = mcc.get_experiments(dataframe=True)
    # get metadata for all non-Cre experiments
    experiments = mca.experiment_source_search(injection_structures = 'root', transgenic_lines = 0)
    ProjPerExp = []
    MetaPerInj = []
    # download the projection density volume for one of the experiments
    for idx,val in enumerate(experiments):
        mca.download_projection_density(infile2, val['id'], resolution = 100)
        MetaPerInj.append({})
        for key,item in val.iteritems():
            MetaPerInj[idx][key] = item
        # read it into memory
        pd_array, pd_info = nrrd.read(infile2)
        ProjPerExp.append(pd_array)
   
    ProjPerExp = np.asarray(ProjPerExp,dtype = 'float32')
    f2 = h5py.File('ProjPerExp.hdf5','w')
    f2.create_dataset('dataset1',data = ProjPerExp)
    pk.dump(MetaPerInj,open('py_files/MetaPerInj.pkl','wb'))
    return MetaPerInj, ProjPerExp


def main():

    print 'Commencing cre-line mission'
    CreMeta, ProjPerCre = GetCreLines()
    print 'Cre-line parsing has been completed'
    print 'Commencing rAAv mission'
    WTMeta, ProjPerExp = ReadConnectivityData()
    print 'rAAv parsing has been completed'

    structureIds,structHash = DownloadStructures() 
    with open('structures.csv', "w") as fp:
        M = []
        for sid in structureIds:
            v = structHash[sid]
            M.append([v['id'], v['acronym'], v['name'],
                      v['parent_structure_id'], v['color_hex_triplet']])
        w = csv.writer(fp, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        w.writerow(['id', 'acronym', 'name',
                    'parent_structure_id', 'color_hex_triplet'])
        for line in M:
            w.writerow(line)

    pk.dump(structureIds,open('structureIds.pkl','wb'))
    pk.dump(structHash,open('structHash.pkl','wb'))

    '''InjPerExp = h5py.File('py_files/InjPerExp.hdf5', 'r')['dataset1']
    ProjPerExp = h5py.File('py_files/ProjPerExp.hdf5', 'r')['dataset1']
    InjPerCre = h5py.File('py_files/InjPerCre.hdf5', 'r')['dataset1']
    CreSet = pk.load(open('py_files/MetaPerCre.pkl', 'rb'))
    NCreSet = pk.load(open('py_files/MetaPerInj.pkl', 'rb'))
    connectivityL = h5py.File('py_files/conL.hdf5', 'r')['dataset1']
    connectivityR = h5py.File('py_files/conR.hdf5', 'r')['dataset1']
    NCreSet2 = pk.load(open('py_files/Cre_clf.pkl', 'rb'))
    layerity = pk.load(open('py_files/target_layerity.pkl', 'rb'))'''


    #WTMeta2 = InjectionClassification(CreMeta,WTMeta)
    #print 'Injection classification is complete'
    #CreMeta2,WTMeta3 = CreFiltering(CreMeta,WTMeta2,'all')
    #print 'Layer filtering is complete'
    unionizes_wt_proj = [QueryAPI(UNIONIZE_CON_FMT % (API_PATH,d['id'])) for d in WTMeta]
    print 'unionizing of wild-type is complete'

    print 'structures have successfully been downloaded'
    connectivityL, connectivityR, connectivityDict = CreateConnectivityMatrix(WTMeta,structureIds,structHash,unionizes_wt_proj)
    fp1 = h5py.File('conL.hdf5','w')
    fp2 = h5py.File('conR.hdf5','w')
    fp1.create_dataset('dataset1',data = connectivityL)
    fp2.create_dataset('dataset1',data = connectivityR)
    print 'Connectivity Matrix stage 1 has been created'


    # Cre-data unionization *************************************************#
    unionizes_cre_proj = [QueryAPI(UNIONIZE_CON_FMT % (API_PATH,d['id'])) for d in CreMeta]
    # pk.dump(unionizes_cre_proj, open('unionizes_cre_proj.pkl', 'wb'))
    print 'unionization is complete'
    cre_pr_L, cre_pr_R = CreateConnectivityMatrix(CreMeta, structureIds, structHash, unionizes_cre_proj)
    fp1 = h5py.File('py_files/cre_pr_L.hdf5', 'w')
    fp2 = h5py.File('py_files/cre_pr_R.hdf5', 'w')
    fp1.create_dataset('dataset1', data = cre_pr_L)
    fp2.create_dataset('dataset1', data = cre_pr_R)
    #************************************************************************#


    # Expression Data Parsing
    ExpMeta = QueryAPI(DATA_EXP_SET_QUERY_URL)
    unionizes_exp = [QueryAPI(UNIONIZE_EXP_FMT % (API_PATH,d['id'])) for d in ExpMeta]
    gene_expression = CreateExpressionMatrix(ExpMeta,structureIds,structHash,unionizes_exp)
    fp1 = h5py.File('py_files/G_Exp.hdf5','w')
    fp1.create_dataset('dataset1',data = gene_expression) 
    pk.dump(ExpMeta, open('GeneMeta.pkl', 'wb'))

    # *********** Step 4: Partition injections based on their cre-line *****#
    cre15      = ['Syt6-Cre_KI148', 'Ntsr1-Cre_GN220', 'Sim1-Cre_KJ18',
                'Efr3a-Cre_NO108', 'Chrna2-Cre_OE25', 'A93-Tg1-Cre',
                'Tlx3-Cre_PL56', 'Rbp4-Cre_KL100', 'Rorb-IRES2-Cre',
                'Scnn1a-Tg3-Cre', 'Nr5a1-Cre', 'Sepw1-Cre_NP39',
                'C57BL/6J', 'Emx1-IRES-Cre', 'Cux2-IRES-Cre']

    creCategories = [cre['transgenic-line'] for cre in CreMeta]
    creCategories = list(set(creCategories))

    Affinity = np.asarray([[0,0,0.001,-5.7125],
                           [-0.001,0,0,5.3625],
                           [0,-0.001,0,5.1625]])

    InjCoo  = []
    InjCoo2 = []
    for idx,injection in enumerate(WTMeta):
      coord = injection['injection-coordinates'];
      coord = np.array([coord[0],coord[1],coord[2],1.0]);
      InjCoo.append( Affinity.dot(coord) )
    InjCoo = np.asarray(InjCoo,dtype = 'float32')      # Wild_type coordinates
    for idx,injection in enumerate(CreMeta):
      coord = injection['injection-coordinates'];
      coord = np.array([coord[0],coord[1],coord[2],1.0]);
      InjCoo2.append( Affinity.dot(coord) )

    CreLineDict = OrderedDict()
    for category in cre15:
        cre_members = np.asarray([idx for idx, val in enumerate(CreMeta)\
        if val['transgenic-line'] == category])
        if len(cre_members) > 0:
            CreLineDict[category] = \
            {'ConMat' : cre_pr_R[:, cre_members],\
            'structure-abbrev' : [CreMeta[val]['structure-abbrev'] for val in cre_members],\
            'layer' : [CreMeta[val]['layer'] for val in cre_members],\
            'cell-type' : [CreMeta[val]['Cell Type'] for val in cre_members],\
            'indices' : cre_members,\
            'id'      : [CreMeta[val]['id'] for val in cre_members],\
            'Coordinates' : InjCoo2[cre_members]}


    CreLineDict['wild_type'] = {'ConMat' : conR, \
                                'structure-abbrev' : \
                                [val['structure-abbrev'] for val in WTMeta],\
                                'layer'    : ['inspecific' for idx in range(len(WTMeta))],\
                                'cell-type': ['inspecific' for idx in range(len(WTMeta))],\
                                'id'       : [val['id'] for val in WTMeta],\
                                'Coordinates' : InjCoo}

    pk.dump(CreLineDict, open('CreLineDict.pkl','wb'))
