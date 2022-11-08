import pandas as pd
import numpy as np
import h5py, json
import app_get_connectivity as GJ

class Preprocess():
    ''' load and process data '''

    def __init__(self,
                 bgdata_path='../BasalGangliaData/data/neurons/striatum/fs',
                 hdf5_output_config_file='output_GJ_correct.hdf5',
                 hdf5_input_config_file=None):
        '''
        Constructor.

        Args:
            bgdata_path (str): Path to fs cell directory
            hdf5_output_config_file (str): file name, including path, to output voltage file
            hdf5_input_config_file (str): file name,including path, to gap junctions and synaptic input
        '''
        self.bgdata_path                = bgdata_path
        self.hdf5_output_config_file    = hdf5_output_config_file
        self.hdf5_input_config_file     = hdf5_input_config_file
        
        self.labels  =   ['gpas', 'epas', 'Ra', 'gnaf', 'qnaf', 'gkdr', \
                            'qkdr', 'gkaf', 'qkaf', 'gkas', 'gbk', 'skir', 'dkir']
        
        self.nparams = len(self.labels)
        
        self.name2directory  = {'FS_0': 'str-fs-e160628_FS2-mMTC180800A-IDB-v20210210',
					            'FS_1': 'str-fs-e161024_FS16-mDR-rat-Mar-13-08-1-536-R-v20210210',
					            'FS_2': 'str-fs-e161205_FS1-mBE104E-v20210209',
					            'FS_3': 'str-fs-e161205_FS1-mMTC180800A-IDB-v20210210' }
        
        self.set_mkey2morph()
        
        self.load_output_config_file()
        self.time               = self.xx['time'][0::2]
        self.ncells             = self.xx['metaData/ID'].shape[0]
        self.depol_threshold    = -40e-3
        self.max_dur_depol      = 20e-3
        
        self.set_id2meta()
	    
        self.check_for_depolBlock() # TODO: get using snudda?
        
        if self.hdf5_input_config_file:
            self.gap_junctions = GJ.Preprocess(hdf5_input_config_file=hdf5_input_config_file)
            self.gap_junctions.load_morphologies()
            self.count_number_of_gj_per_cell()
        else:
            self.gap_junctions = None
	    
        self.get_params()
        
        self.get_voltage()
        
        
    def load_output_config_file(self):
        self.xx = h5py.File(self.hdf5_output_config_file, "r")
    
    def load_json(self, fstr):
        with open(fstr, 'r') as handle:
            return json.load(handle)
    
    
    def check_for_depolBlock(self):
        # basically copy paste on Snudda version. TODO: get using snudda?
        dt 			= self.time[1] - self.time[0]
        n_limit 	= int(np.ceil(self.max_dur_depol / dt))
        count       = 0
        data        = {}
        id2meta     = {}
        for cellid in range(self.ncells):

            volt        = self.xx['neurons'][str(cellid)]['voltage/data'][0]
            ctr         = 0
            id_start    = 0
            
            id2meta[cellid] = { 'family':  str(self.xx['metaData/name'][cellid],'utf8'),
		                        'pkey':     str(self.xx['metaData/parameterKey'][cellid],'utf8'),
		                        'mkey':     str(self.xx['metaData/morphologyKey'][cellid],'utf8')
		                        }

            for idx in range(0, len(self.time)):
                if volt[idx*2] > self.depol_threshold:
                    ctr += 1
                else:
                    if ctr > n_limit/2:
                        count +=1
                        # TODO: remove name, pkey and morph from this dict?
                        if not cellid in data:
                            data[cellid] = {'v': volt[::2],
					                        'start_indx': [id_start],
					                        'end_indx': [idx],
					                        'fs_cell':  str(self.xx['metaData/name'][cellid],'utf8'),
					                        'pkey':     str(self.xx['metaData/parameterKey'][cellid],'utf8'),
					                        'mkey':     str(self.xx['metaData/morphologyKey'][cellid],'utf8')}
                        else:
                            data[cellid]['start_indx'].append(id_start)
                            data[cellid]['end_indx'].append(idx)
                    id_start    = idx
                    ctr         = 0

        self.cells_with_db = data
        self.n_cells_with_DB = count
        self.has_db = [1 if i in self.cells_with_db else 0 for i in range(self.ncells)]
    
    def set_id2meta(self):
    
        id2meta = {}
        for cellid in range(self.ncells):
                
                id2meta[cellid] = { 'family':   str(self.xx['metaData/name'][cellid],'utf8'),
		                            'pkey':     str(self.xx['metaData/parameterKey'][cellid],'utf8'),
		                            'mkey':     str(self.xx['metaData/morphologyKey'][cellid],'utf8')
		                            }
        
        self.id2meta = id2meta
    
    def set_mkey2morph(self):
        self.mkey2morph = {}
        for name in self.name2directory.keys():
            self.mkey2morph[name] = self.load_json(  '{}/{}/{}'.format(  
                                                            self.bgdata_path,  
                                                            self.name2directory[name],
                                                            'morphology/morphology_hash_filename.json'
                                                ))
            
    def count_number_of_gj_per_cell(self):
        self.gj_count = np.zeros(self.ncells)
        for cellid in range(self.ncells):
            count = 0
            if cellid in self.gap_junctions.gj: # cells with no GJ are missing from gap_junctions.gj
                for k,l in self.gap_junctions.gj[cellid].items():
                    count += len(l)
                self.gj_count[cellid] = count
        
        # data for overview histogram
        ipos = np.where(np.array(self.has_db)==1)[0]
        ineg = np.where(np.array(self.has_db)==0)[0]
        self.gj_count_for_hist = {  'pos': self.gj_count[ipos], 
                                    'neg': self.gj_count[ineg],
                                    'rug_text':{'pos': [str(i) for i in ipos], 
                                                'neg': [str(i) for i in ineg]
                                                }
                                    }
    
    def get_params(self, threshold=-0.00):

        res         = np.zeros((self.ncells,self.nparams)) # 
        is_spiking  = np.zeros(self.ncells)
        keys        = []
        
        meta        = {}
        
        if self.hdf5_input_config_file:
            len_morph = np.zeros(self.ncells)
            name_id = np.zeros(self.ncells)

        for cellid in range(self.ncells):
            # get meta data.
            name = str(self.xx['metaData/name'][cellid],'utf8')
            mkey = str(self.xx['metaData/morphologyKey'][cellid],'utf8')
            pkey = str(self.xx['metaData/parameterKey'][cellid],'utf8')
            meta[cellid] = {    'name':name,
                                'mkey':mkey,
                                'pkey':pkey,
                                'morp':'{}/{}/morphology/{}'.format(
                                        self.bgdata_path,
                                        self.name2directory[name],
                                        self.mkey2morph[name][mkey]
                                )}
            # length of morphology
            if self.hdf5_input_config_file:
                len_morph[cellid] = self.gap_junctions.morph_data[name][mkey]['Ltot']
                name_id[cellid] = int(name.split('_')[1])
            # extract channel parameters
            params  = self.load_json('{}/{}/parameters.json'.format(	
                                        self.bgdata_path, 
                                        self.name2directory[meta[cellid]['name']]
                                        ))
            keys.append(meta[cellid]['pkey'])
            pp		= params[meta[cellid]['pkey']]
            if len(pp)<20:
                gpas = pp[2]['value']
                epas = pp[3]['value']
                Ra   = pp[5]['value']
                gnaf = pp[9]['value']
                qnaf = 1.8
                gkdr = pp[10]['value']
                qkdr = pp[11]['value']
                gkaf = pp[12]['value']
                qkaf = 3
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
            
            res[cellid] = [gpas, epas, Ra, gnaf, qnaf, gkdr, qkdr, gkaf, qkaf, gkas, gbk, skir, dkir]

            # check if spikes are crossing 0 potential
            v = self.xx['neurons'][str(cellid)]['voltage/data'][0]
            if max(v) > threshold:
                is_spiking[cellid] = 1

        # jitter params (max 2% of range)
        rmax = np.nanmax(res, axis=0)
        rmin = np.nanmin(res, axis=0)
        r_jitter = res + np.random.uniform(0,0.02*(rmax-rmin), size=res.shape)

        # initiate dataframe
        df = pd.DataFrame(r_jitter, columns=self.labels)

        df['is_spiking'] = is_spiking

        # color of lines based on depol block or no depol block 
        df['has_db']    = self.has_db
        
        # number of GJ/cell
        if self.hdf5_input_config_file:
            df['num_gj'] = self.gj_count
            df['mLtot'] = len_morph
            df['family'] = name_id
            labels = self.labels+['num_gj', 'mLtot', 'family', 'is_spiking', 'has_db']
        else:
            labels = self.labels+['is_spiking', 'has_db']
        
        # insert paramkeys (map pkeys to pkey->index)
        u = np.unique(keys)
        pkey_dict = {}
        for i,k in enumerate(u):
            pkey_dict[k] = i
        key_num = [pkey_dict[k] for k in keys]  
        df.insert(loc=0, column='pkey', value=key_num)
        
        dim = [dict(    range=[  df['pkey'].min(),df['pkey'].max()], 
                        label='pkey', 
                        values=df.loc[:,'pkey'],
                        tickvals=[i for i in range(len(u))],
                        ticktext=u
                        )]
        for i,l in enumerate(labels):
            if l=='family':
                dim.append( dict(    range=[df[l].min(), df[l].max()], 
                                label=l, 
                                values=df.loc[:,l],
                                tickvals=[i for i in range(len(self.gap_junctions.morph_data.keys()))],
                                ticktext=list(self.gap_junctions.morph_data.keys())
                                ))   
            else:
                dim.append( dict(   range=[df[l].min(),df[l].max()], 
                                label=l, 
                                values=df.loc[:,l]
                                ))
        dim.append( dict(   range=[0,self.ncells-1], 
                                label='cellid', 
                                values=[i for i in range(self.ncells)]
                                ))
        self.df = df
        self.dim = dim
        self.meta = meta
    
    def get_voltage(self):
        volt = {}
        for cellid in range(self.ncells):
            volt[cellid] = self.xx['neurons'][str(cellid)]['voltage/data'][0]
        self.voltage = volt
            
    
    def dump_hfpy(self):
        self.xx = None
        
