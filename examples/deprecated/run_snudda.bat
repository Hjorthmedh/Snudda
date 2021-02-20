set SIM_NAME=networks/tinySim

snudda init %SIM_NAME% --size 10 --overwrite

snudda place %SIM_NAME%
snudda detect %SIM_NAME% --volumeID Striatum
snudda prune %SIM_NAME%

cp -a config/input-tinytest-v5.json %SIM_NAME%/input.json

snudda input %SIM_NAME% --input %SIM_NAME%/input.json --time 1.0

snudda simulate %SIM_NAME% --time 0.1
