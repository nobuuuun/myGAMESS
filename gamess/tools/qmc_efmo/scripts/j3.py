#!/usr/bin/env python

from basic import addJ3

f=open("../../../templ_files/J3",'r')
J3_templ_lines=f.readlines()
f.close()

addJ3("h2o.wfs2.xml",J3_templ_lines)
