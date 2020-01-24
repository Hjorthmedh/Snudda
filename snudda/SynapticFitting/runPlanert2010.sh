# Create surrogate data
#python3 Planert2010.py


# Fit to surrogate data
python3 Planert2010part2.py DD < input.txt &> output-DD.txt &
python3 Planert2010part2.py DI < input.txt &> output-DI.txt &
python3 Planert2010part2.py ID < input.txt &> output-ID.txt &
python3 Planert2010part2.py II < input.txt &> output-II.txt &

python3 Planert2010part2.py FD < input.txt &> output-FD.txt &
python3 Planert2010part2.py FI < input.txt &> output-FI.txt &

