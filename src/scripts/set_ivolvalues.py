#! /usr/bin/python3
# Python script for setting interval volume values.

import glob
import os
import subprocess
import sys


# *** Subroutines ***

# Set isovalue in dict scalar_data.
def set_scalar_data(scalar_data, testname, filename, dirname, \
                    isovalueI, isovalueII):
    scalar_data[testname] = dict()
    scalar_data[testname]['filename'] = filename
    scalar_data[testname]['dirname'] = dirname
    scalar_data[testname]['isovalueI'] = isovalueI
    scalar_data[testname]['isovalueII'] = isovalueII


# Add offset to each isovalue.
def add_to_each_isovalue(scalar_data, offsetI, offsetII):

    for tdata in scalar_data:
        scalar_data[tdata]['isovalueI'] += offsetI
        scalar_data[tdata]['isovalueII'] += offsetI

        
# Set integer isovalues.
# scalar_data should be a python dictionary.
def set_isovalues_int(scalar_data):

    # Integer isovalues
    set_scalar_data(scalar_data, 'aneurysm_30_255', 'aneurysm.nhdr', 'volvis', 30, 255)
    set_scalar_data(scalar_data, 'engine_70_140', 'engine.nhdr', 'volvis', 70, 140)
    set_scalar_data(scalar_data, 'frog_45_100', 'Frog.nhdr', 'vollib', 45, 100)
    set_scalar_data(scalar_data, 'fuel_70_150', 'fuel.nhdr', 'volvis', 70, 150)
    set_scalar_data(scalar_data, 'lobster_20_150', 'lobster.nhdr', 'volvis', 20, 150)
    set_scalar_data(scalar_data, 'marschnerlobb_80_120', \
                    'marschnerlobb.nhdr', 'volvis', 80, 120)
    set_scalar_data(scalar_data, 'neghip_70_200', 'neghip.nhdr', 'volvis', 70, 200)
    set_scalar_data(scalar_data, 'silicium_80_160', \
                    'silicium.nhdr', 'volvis', 80, 160)

# Set non-integer isovalues. Increment each isovalue by offset.
# (Typically, offset = 0.5, although any other value is permissible.
def set_isovalues_non_int(scalar_data, offsetI, offsetII):
    set_isovalues_int(scalar_data)
    add_to_each_isovalue(scalar_data, offsetI, offsetII)
