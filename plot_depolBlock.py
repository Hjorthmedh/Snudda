import h5py
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

bgdata_path = '../BasalGangliaData/data/neurons/striatum/fs'
name2directory = {	'FS_0': 'str-fs-e160628_FS2-mMTC180800A-IDB-v20210210',
					'FS_1': 'str-fs-e161024_FS16-mDR-rat-Mar-13-08-1-536-R-v20210210',
					'FS_2': 'str-fs-e161205_FS1-mBE104E-v20210209',
					'FS_3': 'str-fs-e161205_FS1-mMTC180800A-IDB-v20210210' }


def run():
	# FS 100: newer network
	#check_FS100(ncells=250, plot=0)
	
	# first network with depolarization block
	#has_depol_block = print_and_plot_depolBlock()
	has_depol_block = [0, 3, 5, 8, 10, 15, 17, 22, 23, 24, 30, 31, 33, 40, 43, 46, 48]
	check_channel_distribution(has_depol_block)
	


def check_channel_distribution(hdb):
	# hdb is a list of cell ids (that has deplorization block - hdb)	
	import json
	
	xx = h5py.File("output_GJ_correct.hdf5", "r")
	
	labels = ['gpas', 'epas', 'Ra', 'gnaf', 'qnaf', 'gkdr', 'qkdr', 'gkaf', 'qkaf', 'gkas', 'gbk', 'skir', 'dkir']
	n = len(labels)
	
	res = np.zeros((50,n)) # len(hdb)
	
	#for i,cellid in enumerate(hdb):
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
	
	# plot--------------------------
	# param plot
	fig,ax = plt.subplots(1,1)
	for i in range(50):
		if i in hdb:
			ax.plot(range(n), r_01[i], '-o', c='red', alpha=1, zorder=100+i)
		else:
			ax.plot(range(n), r_01[i], '-o', c='grey', alpha=0.3)
	
	ax.set_xticks(range(n))
	ax.set_xticklabels(labels,rotation=90)
	fig.savefig('depol_block_channel_distribution.png', dpi=300)
	
	# box plots
	f2, a2 = plt.subplots(1,n, figsize=(12,3), sharey='all', sharex='all')
	rd = r_01[hdb]
	hdb0 = [i for i in range(50) if i not in hdb]
	r0 = r_01[hdb0]
	for i in range(n):
		a2[i].boxplot([rd[:,i],r0[:,i]])
		a2[i].set_title(labels[i])
		
	a2[0].set_yticks([0,0.5,1])
	for i in range(n):
		a2[i].set_xticks([1,2])
		a2[i].set_xticklabels(['1', '0'])
		a2[i].get_xticklabels()[0].set_color('red')
		a2[i].get_xticklabels()[1].set_color('grey')
	f2.savefig('depol_block_boxplot.png', dpi=300)
    
	bad = {}
	print('------------')
	pdict = {}
	netw = h5py.File("output_GJ.hdf5", "r")
	for i in range(50):
		if r_01[i,2] < 0.35:	
			name 	= str(xx['metaData/name'][i],'utf8')
			fstr 	= '{}/{}/parameters.json'.format(	bgdata_path, 
														name2directory[name])
			with open(fstr, 'r') as handle:
				params = json.load(handle)
			pkey 	= str(xx['metaData/parameterKey'][i],'utf8')
			pp		= params[pkey]
			Ra   	= pp[5]['value']
			plt.figure()
			plt.plot(xx['neurons'][str(i)]['voltage/data'][0])
			plt.title('id={}, Ra01={:.2f}, Ra={:.2f}'.format(i, r_01[i,2], Ra))
			plt.ylim([-.09,.03])
			plt.savefig('Ra_below_150_{}.png'.format(i))
			
			print(name2directory[name], pkey)
			
			if pkey not in pdict:
				for j in range(250):
					pk2 = str(netw['metaData/parameterKey'][j],'utf8')
					if pk2 == pkey:
						plt.figure()
						plt.plot(netw['neurons'][str(j)]['voltage/data'][0])
						plt.title('FS250 cellid:{}'.format(j))
				pdict[pkey] = 0
				plt.show()
	
	print(np.min(res[:,2]), np.max(res[:,2]))		
	plt.show()


def check_FS100(plot=False, ncells=100):
	if ncells == 100:
	    netw = h5py.File("output_100FS.hdf5", "r")
	else:
	    netw = h5py.File("output_GJ.hdf5", "r")    
	time = netw['time'][0:]

	nspikes = np.zeros(ncells)
	ids = [i for i in range(ncells)]

	# loop over cells and check frequency
	for key in range(ncells):
		volt = netw['neurons'][str(key)]['voltage/data'][0]	
		
		if plot:
			plt.plot(time,volt)
			plt.ylim([-0.09,0.03])
			plt.title('cellID={}, nspikes={}'.format(key,len(netw['neurons'][str(key)]['spikes']['data'][0:][0])))
			plt.show()
		
		nspikes[key] = len(netw['neurons'][str(key)]['spikes']['data'][0:][0])
		
	plt.plot(ids, nspikes, 'o', ms=10)
	plt.savefig('spike_freq_pop100FS.png', dpi=300, transparent=True)
	plt.show()
	
	
def print_and_plot_depolBlock(plot_individual = 1):
	xx 			= h5py.File("output_GJ.hdf5", "r") # h5py.File("output_GJ_correct.hdf5", "r")
	time 		= xx['time'][0:]
	dt 			= time[1] - time[0]
	threshold	= -40e-3
	max_dur		= 20e-3
	n_limit 	= int(np.ceil(max_dur / dt))
	fig,ax = plt.subplots(1,1, figsize=(6,2))
	count = 0
	cells_with_depol = []
	name = []
	param = []
	morph = []
	data = {}
	for key in range(250):	#[0,3,8,15,22,24,30,33,46,48]:
		
		volt = xx['neurons'][str(key)]['voltage/data'][0]
		#ax.plot(time, volt)
		
		ctr = 0
		id_start = 0
		
		has_block = []
		
		flag = 0
		for idx in range(0, len(time)):
			if volt[idx] > threshold:
				ctr += 1
			else:
				if ctr > n_limit:
					ax.plot([time[id_start], time[idx]], [count,count], lw=2, c='k')
					if plot_individual:
						if not flag: 
							plt.figure()
							plt.plot(time, volt)
						plt.plot(time[id_start:idx], volt[id_start:idx], '*', color='r')
					flag += 1
				id_start = idx
				ctr = 0

		if flag:
			count += 1
			fs_name = str(xx['metaData/name'][key],'utf8')
			if fs_name not in data:
				data[fs_name] = []
			cells_with_depol.append(key)
			name.append(fs_name)
			param.append(str(xx['metaData/parameterKey'][key],'utf8'))
			morph.append(str(xx['metaData/morphologyKey'][key],'utf8'))
			data[fs_name].append(param[-1])
			if plot_individual:
				plt.ylim([-0.08,0.03])
				plt.title('cell:{}, nblock:{}'.format(key, flag))
			print('ID:{}, name:{}, morph:{}, param:{}'.format( key,	
														xx['metaData/name'][key],
														xx['metaData/morphologyKey'][key],
														xx['metaData/parameterKey'][key] ))
		

	fig.savefig('depol_block_full_population.png', dpi=300, transparent=True)
	#plt.close('all')
	plt.show()	

	print(cells_with_depol)
	print(Counter(name))
	print(Counter(param))
	print(Counter(morph))
	print(data)
	for key in data:
		print('\t{}, {}'.format(key,Counter(data[key])))

	return cells_with_depol



run()
# snudda_load --listN --detailed filen.hdf5 


# get channel profiles of cells ---------------------------------------





