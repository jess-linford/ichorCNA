#!/usr/bin/env python
# coding: utf-8

# In[1]:

import glob
import pathlib
import os
import sys

bamfolder = sys.argv[1]


# In[2]:

bam_files = []
for filepath in pathlib.Path(bamfolder).glob('**/*.bam'):
    bam_files.append(str(filepath.absolute()))
bam_files

# Sort bam files alphabetically by filename
bam_files.sort(key=lambda x: os.path.basename(x))


# In[3]:

f = open("samples.yaml", "w")
f.write('samples:\n')
for bam_file in bam_files:
    f.write(" "+os.path.basename(bam_file)+": "+bam_file+"\n")
f.close()

