# build and run 
python3 setup.py build_ext --inplace

# python3 spkmeans.py 0 wam input_1.txt
# python3 spkmeans.py 0 ddg input_1.txt
# python3 spkmeans.py 2 lnorm input_1.txt
# python3 spkmeans.py 2 jacobi input_1.txt
# python3 spkmeans.py 2 spk input_1.txt

python3 spkmeans.py 0 wam test.csv
python3 spkmeans.py 0 ddg test.csv
python3 spkmeans.py 0 lnorm test.csv
python3 spkmeans.py 0 jacobi Lnorm.txt
python3 spkmeans.py 0 spk test.csv


