#!/bin/python
import argparse

#-------------------------------------------------------------------
# Purpose: Given a list of intent-ins, inouts and outs, this script
# --------
#         print calls for writing variables for the yaml/python I/O
#         generation code
#-------------------------------------------------------------------

def main():

    #parse command line arguments
    ins,outs,finp,fout = parse_args()

    #Break the inputs string into individual variables and remove white spaces
    in_vars = [el.strip() for el in ins.split(',')]
    out_vars = [el.strip() for el in outs.split(',')]



    #Generate write line fortran codes
    #input calls
    write_calls(in_vars, finp)

    #output calls
    write_calls(out_vars, fout, False)

# Rest of the functions

def write_calls(var_list, fname, is_input=True):

    #if the calls are for writing input YAML files,
    # we need both the input and output units as
    #input is also written to the python files
    #Note: commas are the part of these strings
    if is_input:
        unit1 = 'unit_input,'
        unit2 = 'unit_output,'
    #if the calls are for writing output Python files,
    # we need just the output unit as
    #output is only written to the python file
    else:
        unit1 = 'unit_output,'
        unit2 = ''

    #loop through variables and write code for the calls
    with open(fname,'w') as fw:
        for var in var_list:
            fw.write(f"        call write_var({unit1}{unit2}'{var}',{var})\n")

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help="A comma separated list of intent-in/inout variables", required=True)
    parser.add_argument('-o', type=str, help="A comma separated list of intent-out/inout variables", required=True)
    parser.add_argument('--finp', type=str, help="File name for in calls", required=True)
    parser.add_argument('--fout', type=str, help="File name for out calls", required=True)

    args = parser.parse_args()

    #return the args captured by -i and -o flags
    return args.i, args.o, args.finp, args.fout


main()
