#!/bin/bash
#Sets general global environment variable and virtual environment for jobs
echo "Setting up environment"
#Activate virual environement
source ../../../../../snudda_env/bin/activate
####################################################################
#Probably most frequency changed variables
####################################################################
#export NETWORK_NAME=inptun_dspn1
#export NETWORK_NAME=inptun_PD2_dspn1_05 #str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026
export NETWORK_NAME=inptun_PD2_dspn2_01 #str-dspn-e150917_c6_D1-m21-6-DE-v20211028
#export NETWORK_NAME=inptun_PD2_dspn3_01 #str-dspn-e150917_c9_d1-mWT-1215MSN03-v20211026
#export NETWORK_NAME=inptun_PD2_dspn4_01 #str-dspn-e150917_c10_D1-mWT-P270-20-v20211026
#export NETWORK_NAME=inptun_PD2_ispn1_01 #str-ispn-e151123_c1_D2-mWT-P270-09-v20211026
#export NETWORK_NAME=inptun_PD2_ispn2_01 #str-ispn-e150908_c4_D2-m51-5-DE-v20211026
#export NETWORK_NAME=inptun_PD2_ispn3_01 #str-ispn-e160118_c10_D2-m46-3-DE-v20211026
#export NETWORK_NAME=inptun_PD2_ispn4_01 #str-ispn-e150917_c11_D2-mWT-MSN1-v20211026


export NETWORK_DIR=networks/$NETWORK_NAME
export PROJECT_DIR=/cfs/klemming/scratch/${USER:0:1}/$USER/Projects/SnuddaProj09
####################################################################
export SNUDDA_DATA_TOOLS=BasalGangliaData/tools
#SNUDDA_DATA_TOOLS can sometimes be different from $bgdata/tools, so it's better to use separate variable
#bgdata=bgd01/data
bgdata=bgd01/Parkinson/20211105/PD2
if [[ -d "$bgdata" ]]; then
    export SNUDDA_DATA=$bgdata
    echo "Setting SNUDDA_DATA to $SNUDDA_DATA"
else
    echo "SNUDDA_DATA environment variable not changed (may be empty): $SNUDDA_DATA"
fi

####################################################################
#More specific, manual "iterative" variables
#For example if you are tuning cells, you'll want to switch cell folders but don't want to have to do it separately for setup.job, simulate.job, analyse.job etc
####################################################################
export neuronType=dspn
#export neuronType=ispn

#source=bgmod/models/optim/HBP-2021Q4/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026
#source=bgmod/HBP-2021Q4/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026
#source=bgmod/HBP-2021Q4/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026

#############
#source=bgmod/parkinsonian/str-dspn-e150602_c1_D1-mWT-0728MSN01-pd2-v20220114
source=bgmod/parkinsonian/str-dspn-e150917_c6_D1-m21-6-DE-pd2-v20220114
#source=bgmod/parkinsonian/str-dspn-e150917_c9_D1-mWT-1215MSN03-pd2-v20220114
#source=bgmod/parkinsonian/str-dspn-e150917_c10_D1-mWT-P270-20-pd2-v20220114
#source=bgmod/parkinsonian/str-ispn-e151123_c1_D2-mWT-P270-09-pd2-v20220114
#source=bgmod/parkinsonian/str-ispn-e150908_c4_D2-m51-5-DE-pd2-v20220114
#source=bgmod/parkinsonian/str-ispn-e160118_c10_D2-m46-3-DE-pd2-v20220114
#source=bgmod/parkinsonian/str-ispn-e150917_c11_D2-mWT-MSN1-pd2-v20220114
#############

#############
#export singleNeuronType=str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026
export singleNeuronType=str-dspn-e150917_c6_D1-m21-6-DE-v20211028
#export singleNeuronType=str-dspn-e150917_c9_d1-mWT-1215MSN03-v20211026
#export singleNeuronType=str-dspn-e150917_c10_D1-mWT-P270-20-v20211026
#export singleNeuronType=str-ispn-e151123_c1_D2-mWT-P270-09-v20211026
#export singleNeuronType=str-ispn-e150908_c4_D2-m51-5-DE-v20211026
#export singleNeuronType=str-ispn-e160118_c10_D2-m46-3-DE-v20211026
#export singleNeuronType=str-ispn-e150917_c11_D2-mWT-MSN1-v20211026
#############

#############
#export numInputMax=350 #str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026
#export numInputMin=100
export numInputMax=350 #str-dspn-e150917_c6_D1-m21-6-DE-v20211028
export numInputMin=100
#export numInputMax=650 #str-dspn-e150917_c9_d1-mWT-1215MSN03-v20211026
#export numInputMin=400
#export numInputMax=350 #str-dspn-e150917_c10_D1-mWT-P270-20-v20211026
#export numInputMin=100
#export numInputMax=300 #str-ispn-e151123_c1_D2-mWT-P270-09-v20211026
#export numInputMin=50
#export numInputMax=300 #str-ispn-e150908_c4_D2-m51-5-DE-v20211026
#export numInputMin=50
#export numInputMax=260 #str-ispn-e160118_c10_D2-m46-3-DE-v20211026
#export numInputMin=10
#export numInputMax=400 #str-ispn-e150917_c11_D2-mWT-MSN1-v20211026
#export numInputMin=150
#############

export inputType=cortical
#export inputType=thalamic
#export inputType=corticalthalamic #used as dummy to name figures now, but perhaps implement later to tune both types of inputs in one go (perhaps use useMeta)
export useMeta=0

echo "Environment setup at end"
