import getopt
import math
import sys
import argparse
import os
import csv


"""
Util file for parsing text files easily
	
 	- file reads in a .txt file whose separator use space chars: 
 		***specifically PDB and HBondFinder files
	
 	- file may contain a header - if so, use 'header' bool variable and 'header_size' to specify size of header
  
	Usage: parse_file( [file_name] , [optional: header - default = False] , [optional: header_size - default = 1] )

		file_name = name of file to read
		header = boolean indicating if a header is present, and will remove header lines. If header is false, header_size will not be used.
		header_size = number of lines at the top of the file that header occupies.
"""

## Main function to use in this file
def parse_file(file, header=False, header_size=1):
    lines = read_file(file)
    all_lines = parse_lines(lines, header, header_size)
    return all_lines


## Alternate main function specifically for cleaning up PDB files
def parse_PDB_file(file):
    cleaned_lines = read_PDB_lines(file)
    all_lines = parse_lines(cleaned_lines, False, 0)
    return all_lines


def help():
    print("File parses .txt files of data whose separator is empty space characters.")
    print("----------------------------------------------------------------")
    print(
        "Usage: parse_file( [file_name] , [optional: header - default = False] , [optional: header_size - default = 1] )"
    )
    print("")
    print("file_name = name of file to read")
    print(
        "header = boolean indicating if a header is present, and will remove header lines. If header is false, header_size will not be used."
    )
    print("header_size = number of lines at the top of the file that header occupies.")


def read_file(file):
    f = open(file, "r")
    lines = f.readlines()
    return lines


## Utility function to get length of header for PDB files specifically
def get_PDB_header_length(file):
    f = open(file, "r")
    count = 0
    while True:
        line = f.readline()
        if not line:
            break
        if line[0:4] == "ATOM":
            return count
        count += 1


def read_PDB_lines(file):
    f = open(file, "r")
    cleaned_lines = []
    while True:
        line = f.readline()
        if not line:
            break
        if line[0:4] == "ATOM":
            cleaned_lines.append(line)
    return cleaned_lines


def parse_lines(lines, header, header_size):
    all_lines = []
    counter = 0
    if header == False:
        counter = 0
    elif header == True:
        counter = header_size

    for i in range(counter, len(lines)):
        split_line = lines[i].split(" ")
        split_line = [i for i in split_line if i]
        split_line.pop()
        all_lines.append(split_line)
    return all_lines


def main():
    # print(parse_file("hbondfinder_results/HBondFinder_1A4Y.txt"))
    # data = parse_file("Dataset/F_chain_only+h.pdb")
    data = parse_file(
        "Dataset/1brs_barnase_A+h.pdb",
        True,
        get_PDB_header_length("Dataset/1brs_barnase_A+h.pdb"),
    )
    print(len(data))
    data = parse_PDB_file("Dataset/1brs_barnase_A+h.pdb")
    print(len(data))
    # get_PDB_header_length("Dataset/F_chain_only+h.pdb")
    # get_PDB_header_length("Dataset/1brs_barnase_A+h.pdb")


if __name__ == "__main__":
    main()
