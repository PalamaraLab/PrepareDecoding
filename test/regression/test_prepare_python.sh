#!/bin/bash

# Create all decoding quantities files
python3 test_prepare.py

# Test against existing regression output
python3 test_file_contents.py dq_py_original_files.csfs reference_output.csfs || exit 1
python3 test_file_contents.py dq_py_original_files.intervalsInfo reference_output.intervalsInfo || exit 1
python3 test_file_contents.py dq_py_original_files.decodingQuantities.gz reference_output.decodingQuantities.gz || exit 1

echo "Successfully checked regression files"

# Test built-ins with/without frequencies from file
python3 test_file_contents.py dq_py_builtin_except_freq.csfs dq_py_builtin_inc_freq.csfs || exit 1
python3 test_file_contents.py dq_py_builtin_except_freq.intervalsInfo dq_py_builtin_inc_freq.intervalsInfo || exit 1
python3 test_file_contents.py dq_py_builtin_except_freq.decodingQuantities.gz dq_py_builtin_inc_freq.decodingQuantities.gz || exit 1

echo "Successfully checked regression files"

# Tidy up test files
rm dq_py_original_files.csfs
rm dq_py_original_files.intervalsInfo
rm dq_py_original_files.decodingQuantities.gz

rm dq_py_builtin_except_freq.csfs
rm dq_py_builtin_except_freq.intervalsInfo
rm dq_py_builtin_except_freq.decodingQuantities.gz

rm dq_py_builtin_inc_freq.csfs
rm dq_py_builtin_inc_freq.intervalsInfo
rm dq_py_builtin_inc_freq.decodingQuantities.gz
