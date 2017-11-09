# reverse-de-bruijn

ReverseCAKE - Reverse de Bruijn Sequence to Cover All K-mers

ReverseCAKE (short for Reverse de Bruijn sequence to Cover All K-mErs) is a software for generating a shortest sequence that for each k-mer includes the k-mer or its reverse.

Input and Output

ReverseCAKE takes as input four parameters: 
1. k - the order of the sequence. 
2. The alphabets.
3. The output file name. 
4. ILP time limit - if 0, not ILP is run. Need Gurobi library set up to run.

It outputs the sequence as a textual file.

ReverseCAKE was developed by Yaron Orenstein at Ben-Gurion University.

Get the software

Java executable distribution and example files
This distribution is our officially supported executable for ReverseCAKE. This binary is completely self-contained and should work out of the box without any issues.

The software is freely available under the GNU Lesser General Public License, version 3, or any later version at your choice.

ReverseCAKE is a research software, still in the development stage. Hence, it is not presented as error-free, accurate, complete, useful, suitable for any specific application or free from any infringement of any rights. The Software is licensed AS IS, entirely at the user's own risk.

How to use it

java -jar reverse.jar <k> <alphabet> <output_filename> <time limit>

Example runs:

java -jar ReverseCAKE.jar 8 ACGT sequence8.txt 0


Interpreting the output

The output file is a text file containing the sequence.
