#! /usr/bin/python3
# Python script for checking whether ivoldual3 output is a manifold.

import copy
import glob
import os
import subprocess
import sys
from os.path import exists
from set_ivolvalues import set_isovalues_int
from set_ivolvalues import set_isovalues_non_int
from pydos2unix import dos2unix

ijkmeshinfo_prog = 'ijkmeshinfo'
flag_all_tri=True
flag_report_error=True
num_non_manifold = 0
num_self_intersections = 0
flag_integer_isovalues = False
isovalue_offsetI = 0.5
isovalue_offsetII = 0.5
outA_off = "outA.off"
outB_off = "outB.off"
default_test_list = [ 'aneurysm_30_255', 'engine_70_140',\
                      'frog_45_100', 'fuel_70_150', 'lobster_20_150', \
                      'marschnerlobb_80_120',\
                      'neghip_70_200', 'silicium_80_160' ]

# Supersampling allowed only on these data sets.
supersample_test_list = [ 'fuel_70_150', 'neghip_70_200', 'silicium_80_160' ]


# *** Generic subroutines ***

# Run program.
# Note: filepath should be also included in arg_list.
def run_prog(prog_name, arg_list, filepath, flag_report_error):

    global report_file

    cmd = [ prog_name ] + arg_list
    print(' '.join(cmd))
    sys.stdout.flush()
    subprocess_result = subprocess.run(cmd, stdout=report_file)

    if (flag_report_error and (subprocess_result.returncode != 0)):
        print("Error detected in running program " + prog_name + " on file " + filepath + ".")
        print("Exiting...")
        exit(-1)

    return subprocess_result


# Run program.
# Note: filepathA and filepathB should be also included in arg_list.
def run_progII(prog_name, arg_list, filepathA, filepathB, flag_report_error):

    global report_file

    cmd = [ prog_name ] + arg_list
    print(' '.join(cmd))
    sys.stdout.flush()
    subprocess_result = subprocess.run(cmd, stdout=report_file)

    if (flag_report_error and (subprocess_result.returncode != 0)):
        print("Error detected in running program " + prog_name + \
              " on files " + filepathA + " and " + filepathB + ".")
        print("Exiting...")
        exit(-1)

    return subprocess_result


# Convert dos file to unix file
def convert2unix(filename):
    with open(filename, "rb") as src:
        buffer = dos2unix(src)
    with open(filename, "wb") as dest:
        dest.write(buffer)


# *** Run tests Subroutines ***

# Check oriented manifold condition.
def check_manifold(filepath):

    global num_non_manifold
    global num_self_intersections
    
    arg_list = [ "-mesh_dim", "3", "-oriented_manifold", \
                 "-terse", "-report_deep", filepath ]
    result = run_prog(ijkmeshinfo_prog, arg_list, filepath, False)
    if (result.returncode != 0):
        num_non_manifold += 1
    print()
    
    
def check_ivoldual_manifold_on_options \
        (ivoldual_prog, arg_list, isovalueI, isovalueII, filepath):

    global min_diff
    global num_diff
    global flag_check_manifold

    arg_listA = arg_list + \
        [ "-s", "-o", outA_off, str(isovalueI), str(isovalueII), filepath]
    run_prog(ivoldual_prog, arg_listA, filepath, flag_report_error)

    # Check oriented manifold condition.
    check_manifold(outA_off)

    
# *** check ivoldual data structure ***

class check_ivoldual_datastruct:

    def __init__(self, ivoldual_prog, isovalueI, isovalueII, filepath):
        self.prog = ivoldual_prog
        self.isovalueI = isovalueI
        self.isovalueII = isovalueII
        self.filepath = filepath
        self.arg_listA = []

    # Set self.arg_listA
    def SetArgList(self,arg_list):
        self.arg_listA = arg_list

    # Add to self.arg_listA
    def AddToArgList(self,arg_list):
        self.arg_listA += arg_list

    # Check if ivoldual produces an oriented manifold.
    def check_manifold(self):
        check_ivoldual_manifold_on_options\
            (self.prog, self.arg_listA,\
             self.isovalueI, self.isovalueII, self.filepath)
        
    # Check if ivoldual produces an oriented manifold
    #   with options self.arg_listA + arg_listB
    def check_manifoldI(self, arg_listB):
        arg_list = self.arg_listA + arg_listB
        check_ivoldual_manifold_on_options\
            (self.prog, arg_list, \
             self.isovalueI, self.isovalueII, self.filepath)


# Check ivoldual on given test and isovalues.
def check_ivoldual_on_isovalues \
        (ivoldual_prog, ivoldual_options, testname,\
         isovalueI, isovalueII, filename, filepath):

    check_manifold_func = \
        check_ivoldual_datastruct \
        (ivoldual_prog, isovalueI, isovalueII, filepath)

    check_manifold_func.SetArgList(ivoldual_options)
    check_manifold_func.check_manifold()
    check_manifold_func.check_manifoldI(["-trimesh"])

# ivoldual program does not allow call to subsample or supersample with -subdivide.    
#    check_manifold_func.check_manifoldI(["-subsample", "2"])
#    if (testname in supersample_test_list):
#        check_manifold_func.check_manifoldI(["-supersample", "2"])

        
# Check manifold condition on ivoldual run on data.
def check_ivoldual_manifold_on_data\
        (ivoldual_prog, ivoldual_options, testname, testdata):
    
    print()
    print("*** Processing file: ", testdata['filename'], \
          "  Isovalues: ", \
          testdata['isovalueI'], " ", testdata['isovalueII'], " ***")
    print()

    filename = testdata['filename']
    pathname = testdata['dirname'] + os.sep + testdata['filename']
    isovalueI = testdata['isovalueI']
    isovalueII = testdata['isovalueII']

    check_ivoldual_on_isovalues \
        (ivoldual_prog, ivoldual_options, testname, \
         isovalueI, isovalueII, filename, pathname)
                     
    print()


def check_ivoldual_manifold(ivoldual_prog, ivoldual_options, scalar_data, \
                            test_list, report_filename):

    global report_file

    if (report_filename == ""):
        report_file = sys.stdout
        for testname in test_list:
            if testname in scalar_data:
                check_ivoldual_manifold_on_data \
                    (ivoldual_prog, ivoldual_options,
                     testname, scalar_data[testname])
            else:
                print("Error. Unable to find test: ", testname)
                
    else:
        print("Creating report file: " + report_filename)

        with open(report_filename, 'w') as report_file:
            sys.stdout = report_file

            for testname in test_list:
                if testname in scalar_data:
                    check_ivoldual_manifold_on_data \
                        (ivoldual_prog, ivoldual_options, \
                         testname, scalar_data[testname])
                else:
                    print("Error. Unable to find test: ", testname)

            report_summary()

        convert2unix(report_filename)

        sys.stdout = sys.__stdout__


def report_summary():

    if (num_non_manifold > 0):
        print("*** Found", num_non_manifold,
              "instances where output is not an oriented manifold.")
    else:
        print("Passed all manifold checks.")


# *** Usage Subroutines ***
    
# Print usage message.
def usage_msg():
    print("Usage: " + sys.argv[0] + " [OPTIONS] {prog1} [IVOLDUAL_OPTIONS]")
    print("OPTIONS:")
    print("  [-testnames \"list of tests\"]")
    print("  [-integer_isovalues]")
    print("  [-o {report filename}] [-usage]")

    
# Usage error
def usage_error():
    usage_msg();
    exit(-1)

# List testnames
def list_testnames(tlist):
    print("Test names:")
    for testname in tlist:
        print("  ", testname)
        

# *** Main ***

if __name__ == "__main__":

    test_list=default_test_list
    report_test_filename=""

    command_line_argument = list(sys.argv);
    command_line_argument.pop(0);
    
    while ((len(command_line_argument) > 0) and (command_line_argument[0][0] == '-')):
        arg0 = command_line_argument[0]
        if (arg0 == '-testnames'):
            test_list = command_line_argument[1].split()
            command_line_argument.pop(0)
        elif (arg0 == '-trimesh_uniform'):
            flag_trimesh_uniform = True
            flag_all_tri = False
        elif (arg0 == '-no_tri'):
            flag_all_tri = False
            flag_trimesh = False
            flag_trimesh_uniform = False
        elif (arg0 == '-integer_isovalues'):
            flag_integer_isovalues = True
        elif (arg0 == '-o'):
            report_test_filename = command_line_argument[1]
            command_line_argument.pop(0)
        elif (arg0 == '-usage'):
            usage_msg()
            exit(0)
        elif (arg0 == '-list'):
            list_testnames(default_test_list)
            exit(0)
        else:
            print('Illegal option: ' + arg0)
            usage_error();

        command_line_argument.pop(0)

    if (len(command_line_argument) < 1):
        print("Usage error. Missing ivoldual executable name.")
        usage_error()

    ivoldual_prog = command_line_argument[0]
    command_line_argument.pop(0)

    ivoldual_options = command_line_argument

    pathA,fileA = os.path.split(ivoldual_prog)
    if (pathA == ""):
        ivoldual_prog = "./" + ivoldual_prog;

    scalar_data = dict()
    if (flag_integer_isovalues):
        set_isovalues_int(scalar_data)
    else:
        set_isovalues_non_int\
            (scalar_data, isovalue_offsetI, isovalue_offsetII)

    check_ivoldual_manifold(ivoldual_prog, ivoldual_options, scalar_data, \
                            test_list, report_test_filename)

    report_summary()

