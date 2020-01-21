set SIM_NAME=networks/tinySim

snudda init %SIM_NAME% --size 10 --overwrite

REM snudda place %SIM_NAME%
REM snudda detect %SIM_NAME% --volumeID Striatum
REM snudda prune %SIM_NAME%

REM cp -a config/input-tinytest-v5.json %SIM_NAME%/input.json

REM snudda input %SIM_NAME% --input %SIM_NAME%/input.json --time 1.0

REM snudda simulate %SIM_NAME% --time 0.1
