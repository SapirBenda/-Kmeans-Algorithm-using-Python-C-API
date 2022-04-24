# build and run 
python3 setup.py build_ext --inplace

python3 spkmeans.py 0 wam inputs_examples/data_points/input_0.txt
python3 spkmeans.py 0 ddg inputs_examples/data_points/input_1.txt
python3 spkmeans.py 0 lnorm inputs_examples/data_points/input_2.txt
python3 spkmeans.py 0 spk inputs_examples/data_points/input_3.txt
python3 spkmeans.py 0 jacobi inputs_examples/sym_matrix/sym_matrix_input_4.txt




