We needed to update an old network-synapses.hdf5 file, so using ChatGPT we generated the following scripts:

parse_old_file.py --> write 'rename_map.csv'

Then we manually edited 'rename_map.csv'

Now you can run:

convert_old_file.py
add_missing_variables.py

Note that you need to update the path to the file if you want to use this. For some old network.

