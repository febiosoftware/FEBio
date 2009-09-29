#!/bin/bash

svnrev=$(svnversion -n)

sed s/[$]WCREV[$]/$svnrev/ FEBio/svnrev_tpl.txt > FEBio/svnrev.h