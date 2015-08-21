#!/usr/bin/env python
"""
    This is one of the shells to the EbeCollector class. This one
    creates a database using data from subfolders containing multiple
    MUSIC events.
"""

from sys import argv, exit
from os import path

try:
    from_folder = path.abspath(argv[1])
except:
    print("Usage: shell from_folder [sub_folder_pattern] [database_filename]")
    exit()

# get optional parameters
if len(argv)>=3:
    database_filename = argv[2]
else:
    database_filename = "results.db"

# call EbeCollector
from EbeCollector import EbeCollector
EbeCollector().collect_music_results(from_folder, database_filename)
