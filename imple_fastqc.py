#!/usr/bin/python
import os.path
import sys

import subprocess
def imple_fastqc(dir_output,input_file):
#   dir_ouput=sys.argv[0]
 #  input_file=sys.argv[1]	
    cmd = "fastqc -o " + dir_output +" -f fastq " + input_file
    print(cmd)
    os.system(cmd)
	
if __name__=="__main__":
    dir_output=sys.argv[1]
    input_file=sys.argv[2]
    imple_fastqc(dir_output,input_file)
