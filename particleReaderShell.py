#!/usr/bin/env python
"""
    This is one of the shells to the EbeCollector class. This one
    creates a database using data from subfolders containing multiple
    UrQMD events.
"""

from sys import argv, exit
from os import path
from ParticleReader import ParticleReader

try:
    particle_database_filename = argv[1]
except:
    print("Usage: shell to_collect_database_name [output_database_name]")
    exit()

if len(argv)>=3:
    analyzed_database_filename = argv[2]
else:
    analyzed_database_filename = "analyzed_particles.db"

particleDB = ParticleReader(particle_database_filename, 
                            analyzed_database_filename)
particleDB.generateAnalyzedDatabase()
