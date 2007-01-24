#!/usr/bin/python
import sys
import os
print "Working Directory = ",os.getcwd()
for key in os.environ.keys():
    print "%s=:%s"%(key,os.environ[key])

