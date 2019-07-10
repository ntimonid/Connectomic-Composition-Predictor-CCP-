                                         Connectomic Composition Predictor       

The goal of this tool is to facilitate the integration of high-throughput brain related data sources in order to predict and complete a more spatially accurate mouse brain mesoconnectome. Based on the premise that gene expression plays a pivotal role in synapse formation during brain development and homeostasis, the quantity and quality of connectivity information present in gene expression atlases should be assessed. A potential approach is to train supervised machine learning algorithms with gene expression data on the task of learning and predicting  structural connectivity patterns.

For the aforementioned purpose, two brain atlas related data modalities were utilized: 

In situ hybridization (ISH) based gene expression dataset containing the expression energy of 3318 genes for ~1030 spatially distinct mouse brain areas [1].
Anatomical tracing based structural connectivity dataset containing the connectivity strength of projections to ~1030 spatially distinct mouse brain areas. The source locations of projections are derived from a series of multiple tracing experiments that either belong to cre-line categories (transgenic mice) or wild-type (non-transgenic). Currently, the total number of source-injection locations are 1397 [2], [4].
Both data modalities were acquired from the Allen Institute for Brain Science. Documentation can be found at: portal.brain-map.org/

Another Allen Institute based useful resource is the mouse connectivity models toolbox that we have utilized throughout our pipeline and as part of the CCP tool [3]. 

Documentation can be found here: mouse-connectivity-models.readthedocs.io/en/latest/

The supervised machine learning algorithms used to train our models are: Penalized (L2 ridge) linear regression and Random forest regressors.

The CCP tool corresponds to two main classes and a series of independent functions. The two main classes are:

MesoconnectomePredictor: provides preprocessing, storage, analysis and predictive routines for the aforementioned data modalities.
BrainPlotter: provides multiple forms of  visualization for brain related data at different stages of the analysis [5].
Descriptions and implementations of the aforementioned modules can be found in: PrimaryLibrary

One big advantage of the CCP tool is that it provides a user with the capability to visualize their data at any point of their analysis, as well as make new associations between the corresponding data modalities using the predictive and dictionary decomposition routines.

We have provided three independent use cases that indicate ways in which the CCP tool can be utilized for assisting different types of researchers throughout their analysis.

Use case #1: Laminar-resolved Connectomics (see UseCase1.ipynb for source code) shows the ways in which the CCP tool can assist computational neuroscientists.

Use case #2: Predictive Transcriptomics (see UseCase2.ipynb for source code) shows how the CCP tool can assist translational neuroscientists. 

Use case #3 Brain Visualization (see UseCase3.ipynb for source code) shows ways according to which the CCP tool can assist a broader range of people such as high school students, working on projects regarding brain visualization.

The CCP Pipeline (see CCP_Pipeline.ipynb for source code) serves as an archive of all the analytic steps that have been followed in order to develop, test and inspect our neuroinformatics based predictive pipeline.

Documentation regarding description of the pipeline, analysis and interpretation of the results is ongoing.

The Allen_API_Library (see Allen_API_Library.ipynb for source code) is a collection of functions that we have used to download and parse the aforementioned data modalities from the Allen Institute. Due to run-time constraints of our demo, all data were downloaded using this library and were then stored in either .hdf5 or .pkl format.  Throughout the CCP pipeline and all use cases the data were loaded from the storage. Therefore this Library serves as a documentation of how all used data modalities were originally collected and parsed and can be adjusted to the user's needs.

Regarding the execution of jupyter notebook scripts it is important to make a distinction between cells that are on Code or RawNBConvert form. Code implies that the code inside cell is executable while in the second case is just readable. The reason why we have formatted a number of cells as RawNBConvert is because of their corresponding code is quite demanding in terms of execution time and therefore we keep the cell in readable form in order to show the user how the code flows, while the cell below that one serves as a short-circuit solution where the results are being loaded from pickle files where they had been stored rather than re-run. 

References

[1] E.S. Lein, M.J. Hawrylycz, A. Nancy et al. Genome-wide atlas of gene expression in the adult mouse brain. Nature, 445:168-176, 2007.

[2] S.W. Oh, J.A. Harris, L. Ng et al. A mesoscale connectome of the mouse brain. Nature, 508:207-213, 2014.

[3] J.E. Knox et al. High-resolution data-driven model of the mouse connectome. Network Neuroscience, 3 (1): 217-236, 2019.

[4] J.A. Harris et al. The organization of intracortical connections by layer and cell class in the mouse brain.  Biorxiv, 2018.

[5] R. Bakker, P. Tiesinga, R. Kötter. The Scalable Brain Atlas: Instant Web-Based Access to Public Brain Atlases and Related Content. Neuroinformatics, 13 (3): 353–366, 2015.

[6] S. Ji, A. Fakhry, H. Deng. Integrative analysis of the connectivity and gene expression atlases in the mouse brain. Neuroimage 84:245–253  (2014). doi:10.1016/j.neuroimage.2013.08.049.
