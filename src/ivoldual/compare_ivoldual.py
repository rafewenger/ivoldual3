#! /usr/bin/python3
# Python script for comparing ivoldual executables.

import copy
import glob
import os
import subprocess
import sys
from os.path import exists
from set_ivolvalues import set_isovalues_int
from set_ivolvalues import set_isovalues_non_int
from pydos2unix import dos2unix

ijkmeshdiff_prog = 'ijkmeshdiff'
ijkmeshinfo_prog = 'ijkmeshinfo'
flag_all_tri=True
flag_check_manifold=True
flag_report_error=True
flag_sort=False
min_diff = 0.0
num_diff = 0
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
    
    
def compare_ivoldual_on_options(ivoldual_progA, ivoldual_progB, arg_list, \
                                isovalueI, isovalueII, filepath):

    global min_diff
    global num_diff
    global flag_sort
    global flag_check_manifold

    arg_listA = arg_list + \
        [ "-s", "-o", outA_off, str(isovalueI), str(isovalueII), filepath]
    run_prog(ivoldual_progA, arg_listA, filepath, flag_report_error)

    arg_listB = arg_list + \
        [ "-s", "-o", outB_off, str(isovalueI), str(isovalueII), filepath]
    run_prog(ivoldual_progB, arg_listB, filepath, flag_report_error)

    arg_list_mdiff = [ "-terse", outA_off, outB_off ]
    if (min_diff > 0.0):
        arg_list_mdiff = \
            [ "-min_diff", str(min_diff) ] + arg_list_mdiff
    if (flag_sort):
        arg_list_mdiff = [ "-sort" ] + arg_list_mdiff
    
    result = run_progII(ijkmeshdiff_prog, arg_list_mdiff, \
                        outA_off, outB_off, False)
    print()
    
    if (result.returncode != 0):
        num_diff += 1

        if flag_check_manifold:
            # Check oriented manifold condition.
            check_manifold(outA_off)

    
# *** compare ivoldual data structure ***

class compare_ivoldual_datastruct:

    def __init__(self, ivoldual_progA, ivoldual_progB,
                 isovalueI, isovalueII, filepath):
        self.progA = ivoldual_progA
        self.progB = ivoldual_progB
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
        
    # Compare on options self.arg_listA
    def compare(self):
        compare_ivoldual_on_options\
            (self.progA, self.progB, self.arg_listA,\
             self.isovalueI, self.isovalueII, self.filepath)
        
    # Compare on options self.arg_listA + arg_listB
    def compareI(self, arg_listB):
        arg_list = self.arg_listA + arg_listB
        compare_ivoldual_on_options\
            (self.progA, self.progB, arg_list, \
             self.isovalueI, self.isovalueII, self.filepath)

# Compare ivoldual with option -trimesh.
def compare_ivoldual_trimesh(compare_func):

    compare_func.compareI(["-trimesh"])

    
# Compare ivoldual with option -trimesh_uniform.
def compare_ivoldual_trimesh_uniform(compare_func):
    compare_func.compareI(["-trimesh_uniform"])
    
        
# Compare ivoldual.
def compare_ivoldual(ivoldual_progA, ivoldual_progB, testname,\
                     isovalueI, isovalueII, filename, filepath):

    compare_func = \
        compare_ivoldual_datastruct\
        (ivoldual_progA, ivoldual_progB, isovalueI, isovalueII, filepath)

    compare_func.compare()
    compare_func.compareI(["-trimesh"])
    compare_func.compareI(["-subdivide"])
    compare_func.compareI(["-subdivide", "-trimesh"])
    compare_func.compareI(["-split_hex"])
    compare_func.compareI(["-subsample", "2"])
    if (testname in supersample_test_list):
        compare_func.compareI(["-supersample", "2"])

    if (flag_trimesh):
        compare_ivoldual_trimesh(compare_func)

    if (flag_trimesh_uniform):
        compare_ivoldual_trimesh_uniform(compare_func)

        
# Compare programs on data.
def compare_prog_on_data(ivoldual_progA, ivoldual_progB, testname, testdata):
    
    print()
    print("*** Processing file: ", testdata['filename'], \
          "  Isovalues: ", \
          testdata['isovalueI'], " ", testdata['isovalueII'], " ***")
    print()

    filename = testdata['filename']
    pathname = testdata['dirname'] + os.sep + testdata['filename']
    isovalueI = testdata['isovalueI']
    isovalueII = testdata['isovalueII']

    compare_ivoldual(ivoldual_progA, ivoldual_progB, testname,\
                     isovalueI, isovalueII, filename, pathname)
                     
    print()


def compare_prog(ivoldual_progA, ivoldual_progB, \
                 scalar_data, test_list, report_filename):

    global report_file

    if (report_filename == ""):
        report_file = sys.stdout
        for testname in test_list:
            if testname in scalar_data:
                compare_prog_on_data(ivoldual_progA, ivoldual_progB,\
                                     testname, scalar_data[testname])
            else:
                print("Error. Unable to find test: ", testname)
                
    else:
        print("Creating report file: " + report_filename)

        with open(report_filename, 'w') as report_file:
            sys.stdout = report_file

            for testname in test_list:
                if testname in scalar_data:
                    compare_prog_on_data(ivoldual_progA, ivoldual_progB,\
                                         testname, scalar_data[testname])
                else:
                    print("Error. Unable to find test: ", testname)

            report_summary()

        convert2unix(report_filename)

        sys.stdout = sys.__stdout__


def report_summary():

    if (num_diff > 0):
        print("*** Found", num_diff,
              "instances where output files do not match.")

        if (num_non_manifold > 0):
            print("*** Found", num_non_manifold,
                  "instances where output is not an oriented manifold.")
        else:
            print("Passed manifold checks wherever differences were found.");

        if (num_self_intersections > 0):
            print("*** Found", num_self_intersections,
                  "instances where output has self intersections.")
        else:
            print("Passed self intersection checks wherever differences were found.");

    else:
        print("No differences between any output files.");


# *** Usage Subroutines ***
    
# Print usage message.
def usage_msg():
    print("Usage: " + sys.argv[0] + " [OPTIONS] {prog1} {prog2}")
    print("OPTIONS:")
    print("  [-min_diff {min difference}]")
    print("  [-testnames \"list of tests\"]")
    print("  [-trimesh_uniform | -no_tri] [-skip_manifold_check]")
    print("  [-sort]")
    print("  [-integer_isovalues]")
    print("  [-o {report filename}] [-usage]")

    
# Usage error
def usage_error():
    usage_msg();
    exit(-1)


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
        elif (arg0 == '-min_diff'):
            min_diff = float(command_line_argument[1])
            command_line_argument.pop(0)
        elif (arg0 == '-sort'):
            flag_sort = True
        elif (arg0 == '-trimesh_uniform'):
            flag_trimesh_uniform = True
            flag_all_tri = False
        elif (arg0 == '-no_tri'):
            flag_all_tri = False
            flag_trimesh = False
            flag_trimesh_uniform = False
        elif (arg0 == '-skip_manifold_check'):
            flag_check_manifold = False
        elif (arg0 == '-integer_isovalues'):
            flag_integer_isovalues = True
        elif (arg0 == '-o'):
            report_test_filename = command_line_argument[1]
            command_line_argument.pop(0)
        elif (arg0 == '-usage'):
            usage_msg()
            exit(0)
        else:
            print('Illegal option: ' + arg0)
            usage_error();

        command_line_argument.pop(0)

    if (len(command_line_argument) < 1):
        print("Usage error. Missing ivoldual executable names.")
        usage_error()

    if (len(command_line_argument) < 2):
        print("Usage error. Missing second ivoldual executable name.")
        usage_error()

    if (len(command_line_argument) > 2):
        print("Usage error. Expecting exactly two executables names.")
        usage_error()

    ivoldual_progA = command_line_argument[0]
    ivoldual_progB = command_line_argument[1]

    pathA,fileA = os.path.split(ivoldual_progA)
    if (pathA == ""):
        ivoldual_progA = "./" + ivoldual_progA;

    pathB,fileB = os.path.split(ivoldual_progB)
    if (pathB == ""):
        ivoldual_progB = "./" + ivoldual_progB;

    scalar_data = dict()
    if (flag_integer_isovalues):
        set_isovalues_int(scalar_data)
    else:
        set_isovalues_non_int\
            (scalar_data, isovalue_offsetI, isovalue_offsetII)
    
    if (flag_all_tri):
        flag_trimesh = True
        flag_trimesh_uniform = True

    compare_prog(ivoldual_progA, ivoldual_progB, scalar_data, \
                 test_list, report_test_filename)

    report_summary()

