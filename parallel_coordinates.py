import plotly.express as px

import pandas as pd
import plotly.graph_objs as go
import plotly.tools as tls

import numpy as np
import h5py, json

bgdata_path = '../BasalGangliaData/data/neurons/striatum/fs'
name2directory = {	'FS_0': 'str-fs-e160628_FS2-mMTC180800A-IDB-v20210210',
					'FS_1': 'str-fs-e161024_FS16-mDR-rat-Mar-13-08-1-536-R-v20210210',
					'FS_2': 'str-fs-e161205_FS1-mBE104E-v20210209',
					'FS_3': 'str-fs-e161205_FS1-mMTC180800A-IDB-v20210210' }



# pre-process --------------------------------------------------------
xx = h5py.File("output_GJ_correct.hdf5", "r")

labels = ['gpas', 'epas', 'Ra', 'gnaf', 'qnaf', 'gkdr', 'qkdr', 'gkaf', 'qkaf', 'gkas', 'gbk', 'skir', 'dkir']
n = len(labels)

res = np.zeros((50,n)) # len(hdb)

for i,cellid in enumerate([j for j in range(50)]):
	# extract channel parameters
	name 	= str(xx['metaData/name'][cellid],'utf8')
	fstr 	= '{}/{}/parameters.json'.format(	bgdata_path, 
												name2directory[name])
	with open(fstr, 'r') as handle:
		params = json.load(handle)
	pkey 	= str(xx['metaData/parameterKey'][cellid],'utf8')
	pp		= params[pkey]
	if len(pp)<20:
		gpas = pp[2]['value']
		epas = pp[3]['value']
		Ra   = pp[5]['value']
		gnaf = pp[9]['value']
		qnaf = None
		gkdr = pp[10]['value']
		qkdr = pp[11]['value']
		gkaf = pp[12]['value']
		qkaf = None
		gkas = pp[13]['value']
		gbk  = pp[15]['value']
		skir = pp[16]['value']
		dkir = pp[17]['value']
	else:
		gpas = pp[2]['value']
		epas = pp[3]['value']
		Ra   = pp[5]['value']
		gnaf = pp[9]['value']
		qnaf = pp[10]['value']
		gkdr = pp[11]['value']
		qkdr = pp[12]['value']
		gkaf = pp[13]['value']
		qkaf = pp[14]['value']
		gkas = pp[15]['value']
		gbk  = pp[17]['value']
		skir = pp[18]['value']
		dkir = pp[19]['value']
	
	res[i] = [gpas, epas, Ra, gnaf, qnaf, gkdr, qkdr, gkaf, qkaf, gkas, gbk, skir, dkir]

# normalize ranges to 0-1
rmax = np.nanmax(res, axis=0)
rmin = np.nanmin(res, axis=0)
r_01 = (res-rmin)/(rmax-rmin)

r_jitter = r_01 + np.random.uniform(0,0.02, size=r_01.shape)

# initiate dataframe
df = pd.DataFrame(r_jitter, columns=labels)

# color of lines based on depol block or no depol block
has_depol_block = [0, 3, 5, 8, 10, 15, 17, 22, 23, 24, 30, 31, 33, 40, 43, 46, 48]
hdb             = [1 if i in has_depol_block else 0 for i in range(50)]
df['hdb']       = hdb

ldict = {}
for l in labels:
    ldict[l] = l


fig = px.parallel_coordinates(df, color="hdb", labels=ldict,
                             color_continuous_scale=px.colors.diverging.Tealrose,
                             color_continuous_midpoint=0.5)


fig.show()

