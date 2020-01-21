simName=networks/tinySim

#./snudda.py init $simName --size 1760000
./snudda.py init $simName --size 10 --overwrite #00

./snudda.py place $simName
#./snudda.py detect $simName --hvsize 50 --volumeID Striatum
./snudda.py detect $simName --volumeID Striatum
./snudda.py prune $simName

# Copy over template input
cp -a config/input-tinytest-v5.json $simName/input.json

# Uncomment this to generate input
 ./snudda.py input $simName --input $simName/input.json --time 1.0

# Uncomment this to run simulation
./snudda.py simulate $simName --time 0.1
