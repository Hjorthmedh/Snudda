import pandas as pd
import numpy as np
import h5py, json
import morph_lib_creator

class Preprocess():
    ''' load and process data '''

    def __init__(self,
                 bgdata_path='../BasalGangliaData/data/neurons/striatum/fs',
                 hdf5_input_config_file='network-synapses.hdf5'):
        '''
        Constructor.

        Args:
            bgdata_path (str): Path to fs cell directory
            hdf5_input_config_file (str): file name, including path, to gap junction file
        '''
        self.bgdata_path                = bgdata_path
        self.hdf5_input_config_file     = hdf5_input_config_file
        
        self.name2directory  = {'FS_0': 'str-fs-e160628_FS2-mMTC180800A-IDB-v20210210',
					            'FS_1': 'str-fs-e161024_FS16-mDR-rat-Mar-13-08-1-536-R-v20210210',
					            'FS_2': 'str-fs-e161205_FS1-mBE104E-v20210209',
					            'FS_3': 'str-fs-e161205_FS1-mMTC180800A-IDB-v20210210' }
        
        self.load_output_config_file()
        
        self.get_gapjunctions()
        
        #self.load_morphologies()
        
        
    def load_output_config_file(self):
        self.xx = h5py.File(self.hdf5_input_config_file, "r")
    
    def get_gapjunctions(self):
        ''' 
        the "gapJunctions" dataset is an array where each row is a gapjunction.
        The position in each row give:
            0: sourceCellID, 1: destCellID, 2: sourceSecID, 3: destSecID,
            4: sourceSegX, 5: destSegX, 6: voxelX, 7: voxelY, 8: voxelZ,
            9: hyperVoxelID, 10: conductance (integer, in pS)
        '''
        # get the dataset. 
        dset = self.xx['network']['gapJunctions'][:,:]
        
        # collect all GJ in single cells by reformating into dict: cellid->secid = [segx1,...,segxN] 
        # (I have controlled so that a pre-post cell-pair is not part of the list as post-pre)
        gj = {}
        for row in dset:    # we do pre-post first
            if row[0] not in gj:
                gj[row[0]] = {row[2]:[row[4]]}
            elif row[2] not in gj[row[0]]:
                gj[row[0]][row[2]] = [row[4]]
            else:
                gj[row[0]][row[2]].append(row[4])
        for row in dset:    # post-pre (perhaps this could be done in the same loop?)
            if row[1] not in gj:
                gj[row[1]] = {row[3]:[row[5]]}
            elif row[3] not in gj[row[1]]:
                gj[row[1]][row[3]] = [row[5]]
            else:
                gj[row[1]][row[3]].append(row[5])
        
        self.gj = gj
        
        
    def load_morphologies(self):
        morph_data = {}
        for key in self.name2directory.keys():
            path = '{}/{}/morphology'.format(   self.bgdata_path, 
                                                self.name2directory[key])
            with open('{}/morphology_hash_filename.json'.format(path), 'r') as handle:
                mkey2fname = json.load(handle)
            morph_data[key] = {}
            for mkey,fname in mkey2fname.items():
                fname = '{}/{}'.format(   path, 
                                          fname)
                morph_with_none, sec_coordinates, stem2plot, sec2stem, morphology = morph_lib_creator.create(
                    swc_file=fname)
                morph_data[key][mkey] = {'morph_with_none':morph_with_none,
                                        'sec_coordinates':sec_coordinates,
                                        'stem2plot':stem2plot,
                                        'sec2stem':sec2stem,
                                        'morphology':morphology,
                                        'Ltot':morph_lib_creator.calc_total_length(morphology['points'])
                                        }
            
        self.morph_data = morph_data
    
    def set_morphology_Ltot(self, points):
        # calculates the total dendritic length of morphologies
        #   morphologies must be loaded first...
        if not morph_data in self:
            self.load_morphologies()
        for cell in self.morph_data:
            for mkey in self.morph_data[cell]:
                p = self.morph_data[cell][mkey]['morphology']['points'] 
                self.morph_data[cell][mkey]['Ltot'] = morph_lib_creator.calc_total_length(p)
    
    
    
    
