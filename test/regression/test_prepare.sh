#!/bin/bash

prepareDecodingExe=$1

if [ -f "$prepareDecodingExe" ]; then
    echo "Using $prepareDecodingExe as prepare decoding executable."
else
    echo "Please specify the prepare decoding executable as argument to this script." && exit 1
fi

# demographics file
demoFile=input_CEU.demo

# Time discretization file
discFile=input_30-100-2000.disc

# number of CSFS individuals. Values >= 50 should be sufficient, 300 is good. Needs to be <= number of (haploid) samples in the data set
CSFSsamples=50

# root of file to be analyzed
outFile=test_output

# frequency file
freqFile=input_UKBB.frq

# precompute decoding quantities
echo "Precomputing decoding quantities"
./"$prepareDecodingExe" \
	-D ${demoFile} \
	-d ${discFile} \
	-n ${CSFSsamples} \
	-F ${freqFile} \
	-o ${outFile} || exit 1
echo "Finished precomputing decoding quantities"

python3 test_file_contents.py test_output.csfs reference_output.csfs || exit 1
python3 test_file_contents.py test_output.intervalsInfo reference_output.intervalsInfo || exit 1
python3 test_file_contents.py test_output.decodingQuantities.gz reference_output.decodingQuantities.gz || exit 1

echo "Successfully checked files"

rm test_output.csfs
rm test_output.intervalsInfo
rm test_output.decodingQuantities.gz
