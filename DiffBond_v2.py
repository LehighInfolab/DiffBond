import getopt
import math
import sys
import argparse
from tempfile import tempdir

import PDBGreedySearch
import PDB_HB_parser

"""
GreedySearch main function
- Takes 2 PDB files, optional distance variable, optional output file NameError
"""

# Parse arguments from command line


def parseArg():
    parser = argparse.ArgumentParser(
        description="Identify all points between protein structures or chains that are within a certain distance of each other."
    )
    parser.add_argument(
        "-i",
        nargs="+",
        metavar="InputPDB",
        help="Input PDB file to be compared. If only 1 file as input, then DiffBond will find bonds between all protein chains. If 2 files as input, then DiffBond will find bonds between the 2 PDB files.",
    )
    parser.add_argument("-o", nargs="?", metavar="OutputPDB", help="Output file name")

    parser.add_argument(
        "-m",
        nargs="+",
        metavar="mode",
        help="Search mode can be multiple combinations of the following options. Must include at least 1 option. Contact = c, Ionic bond = i, Hydrogen bond = h, Salt bridge = S, Cation pi = p",
    )
    parser.add_argument(
        "-d",
        nargs="?",
        metavar="distance",
        type=float,
        help="Resolution for distance checking.",
    )

    # parse list of args
    args = parser.parse_args()
    args = vars(args)
    i_list = args["i"]
    m_list = args["m"]
    d = args["d"]
    o = args["o"]
    return i_list, m_list, d, o


# Checks a connected graph from all atoms in point1 to point2, where both point1 and point2 are lists of parsed values from a PDB file.
# If within dist, then append points1, points2, and dist as a list to output. This function uses the 3D distance equation to calculate distance.
def compareDist(points1, points2, dist):
    output = []
    for i in points1:
        for j in points2:
            d1 = (float(i[6]) - float(j[6])) ** 2
            d2 = (float(i[7]) - float(j[7])) ** 2
            d3 = (float(i[8]) - float(j[8])) ** 2
            d = math.sqrt(d1 + d2 + d3)
            if d < dist:
                edge = [i, j, d]
                output.append(edge)
    return output


## Util function to check the charge of a point based on the residue and atom. This is used for getting ionic bonds.
def get_charge_from_res(res, atom):
    charge = 0
    if res == "HIS" or res == "ARG" or res == "LYS":
        if "NZ" in atom or "NE" in atom or "ND" in atom or "NH" in atom:
            charge = 1
    elif res == "ASP" or res == "GLU":
        if "OD" in atom or "OE" in atom:
            charge = -1
    return charge


# Function for finding atoms meeting ionic bond criteria. Checks a connected graph from charged atoms in point1 to point2.
# If within dist, then append points1, points2, and dist as a list to output. This function uses the 3D distance equation to calculate distance.
def compareDistIonic(points1, points2, dist):
    output = []

    # Swap lists of points1 and points2 if points1 is larger. Slight improvement in efficiency if points1 much larger since function checks point1 charge before going into 2nd for loop.
    if len(points1) > len(points2):
        temp = points1
        points1 = points2
        points2 = temp

    # Temp points lists to hold all points that meet charge criteria before calculating distance.
    p1 = []
    p2 = []
    for i in points1:
        res1 = i[3]
        atom1 = i[2]
        charge1 = get_charge_from_res(res1, atom1)
        if charge1 != 0:

            for j in points2:
                res2 = j[3]
                atom2 = j[2]
                charge2 = get_charge_from_res(res2, atom2)
                if charge2 != 0:

                    # If charge1 and charge2 are not 0, then they must be 1 or -1. If charge1 + charge2 = 0, then one must be -1 and other must be +1.
                    total_charge = charge1 + charge2
                    if total_charge == 0:
                        if i not in p1:
                            p1.append(i)
                        if j not in p2:
                            p2.append(j)

    output = compareDist(p1, p2, dist)
    return output


# TODO: still exploring whether cation pi calculations are possible.
# def compareDistCatPi(output, points1, points2, dist, splitLines_PDB1, splitLines_PDB2):
#     output1 = []
#     for i in range(len(output)):
#         pts1_idx = int(output[i][0])
#         temp_AminoAcid = splitLines_PDB1[pts1_idx][3]
#         charge1 = 0
#         temp = [pts1_idx]
#         temp.append([])
#         if temp_AminoAcid == "ARG" or temp_AminoAcid == "LYS":
#             charge1 = 1
#         elif (
#             temp_AminoAcid == "PHE"
#             or temp_AminoAcid == "TYR"
#             or temp_AminoAcid == "TRP"
#         ):
#             charge1 = -1
#         if charge1 != 0:
#             for j in range(len(output[i][1])):
#                 charge2 = 0
#                 pts2_idx = int(output[i][1][j])
#                 temp_AminoAcid = splitLines_PDB2[pts2_idx][3]
#                 if temp_AminoAcid == "ARG" or temp_AminoAcid == "LYS":
#                     charge2 = 1
#                 elif (
#                     temp_AminoAcid == "PHE"
#                     or temp_AminoAcid == "TYR"
#                     or temp_AminoAcid == "TRP"
#                 ):
#                     charge2 = -1
#                 total_charge = charge1 + charge2
#                 if total_charge == 0:
#                     # d1 = (points1[pts1_idx][0]-points2[pts2_idx][0])**2
#                     # d2 = (points1[pts1_idx][1]-points2[pts2_idx][1])**2
#                     # d3 = (points1[pts1_idx][2]-points2[pts2_idx][2])**2
#                     # d	=	math.sqrt(d1+d2+d3)
#                     # if d < dist:
#                     temp[1].append(str(pts2_idx))
#         if temp[1]:
#             output1.append(temp)
#     return output1


# Util function for Removing duplicates and also parses out any atoms with # = -1 or 0
# Since changing the code to use only line indexes, duplicates will not occur so this function is unused.
def removeDupe(output):
    output = list(dict.fromkeys(output))
    for i in output:
        if i == -1:
            output.remove(-1)
        if i == 0:
            output.remove(0)
    return output


def switch(mode):
    if mode == "c":
        contact_func = lambda a, b, c: compareDist(a, b, c)
        return contact_func
    elif mode == "i":
        ionic_func = lambda a, b, c: compareDistIonic(a, b, c)
        return ionic_func


def main():
    # i_list, mode, dist, outputFileName = parseArg()
    # i_list = ["Dataset/1brs_barnase_A+h.pdb", "Dataset/F_chain_only+h.pdb"]
    # i_list = ["Dataset/F_chain_only+h.pdb", "Dataset/6cnk-F_chain+h.pdb"]
    i_list = ["Dataset/6cnk-F_chain+h.pdb", "Dataset/F_chain_only+h.pdb"]
    # i_list = ["Dataset/F_chain_only+h.pdb", "Dataset/F_chain_only+h.pdb"]
    mode = ["i"]
    dist = None
    outputFileName = None

    if dist == None:
        dist = 5

    if outputFileName == None:
        outputFileName = "Test"
        for i in i_list:
            outputFileName = outputFileName + "_" + str(i.split(".")[0])
        outputFileName = outputFileName + "_" + str(mode)
        outputFileName = outputFileName + "_" + str(dist)

    if len(i_list) == 1:
        print("1")
    elif len(i_list) == 2:
        PDB_data = []
        for i in i_list:
            data = PDB_HB_parser.parse_PDB_file(i)
            PDB_data.append(data)

        for m in mode:
            func = switch(m)
            edges = func(PDB_data[0], PDB_data[1], dist)
            print(edges)


if __name__ == "__main__":
    main()
