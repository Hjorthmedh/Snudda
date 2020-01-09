export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

ipcluster start -n 12 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

# simName=networks/article100000-v1
#simName=networks/res5test
# simName=networks/article2019-v6
#simName=networks/article2019-v13
#simName=networks/SPN2SPN-distDepPruning-v5
#simName=networks/LTS-con-v7
#simName=networks/Robert-bug-v1
#simName=networks/NewConfig-v4
simName=networks/FullStriatum-v2

if [ -d "$simName" ]; then
  echo "Directory $simName already exists!!"
  exit 666    
fi

./snudda.py init $simName --size 1760000
#./snudda.py init $simName --size 100000 #3000 #200000

./snudda.py place $simName 
#./snudda.py detect $simName --hvsize 50 --volumeID Striatum
./snudda.py detect $simName --volumeID Striatum
./snudda.py prune $simName

# Copy over template input
#cp -a config/input-tinytest-v4.json $simName/input.json

# Uncomment this to generate input
# ./snudda.py input $simName --input $simName/input.json

ipcluster stop


# Uncomment this to run simulation
# ./snudda.py simulate $simName

