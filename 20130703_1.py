#!/usr/bin/env python
#!vim:fileencoding=UTF-8


import subprocess

jobid = (
("sf_0002", "A_onlyAICG"),
("sf_0004", "A_onlyAICG"),
("sf_0009", "I_ELE_HIS0_P1all"),
("sf_0010", "I_ELE_HIS0_P1all"),
("sf_0011", "G_ELE_HIS0_noP"),
("sf_0012", "G_ELE_HIS0_noP"),
("sf_0015", "J_ELE_HIS0_P2act"),
("sf_0016", "J_ELE_HIS0_P2act"),
("sf_0017", "K_ELE_HIS0_P2all"),
("sf_0018", "K_ELE_HIS0_P2all"),
("sf_0020", "A_onlyAICG"),
("sf_0021", "A_onlyAICG"),
("sf_0022", "A_onlyAICG"),
("sf_0023", "G_ELE_HIS0_noP"),
("sf_0024", "G_ELE_HIS0_noP"),
("sf_0025", "G_ELE_HIS0_noP"),
("sf_0026", "K_ELE_HIS0_P2all"),
("sf_0027", "K_ELE_HIS0_P2all"),
("sf_0028", "K_ELE_HIS0_P2all"),
("sf_0029", "J_ELE_HIS0_P2act"),
("sf_0030", "J_ELE_HIS0_P2act"),
("sf_0031", "J_ELE_HIS0_P2act"),
("sf_0032", "I_ELE_HIS0_P1all"),
("sf_0033", "I_ELE_HIS0_P1all"),
("sf_0034", "I_ELE_HIS0_P1all"),
("sf_0035", "L"), ("sf_0036", "L"),
("sf_0037", "L"), ("sf_0038", "L"), ("sf_0039", "L"),
("sf_0040", "T"), ("sf_0041", "T"), ("sf_0042", "T"),
("sf_0043", "T"), ("sf_0044", "T"), ("sf_0045", "S"),
("sf_0046", "S"), ("sf_0047", "S"),
)

pathroot = "/home/hori/mapk/cafemol/"

for job in jobid:
    jobname = job[0]
    group = job[1]
    
    wd = pathroot + jobname
    cmdline = "20130702_3.py polar_f3.out " + jobname
    p = subprocess.Popen(cmdline, shell=True, cwd=wd)
    p.wait()
    
    cmdline = "gnuplot ../hist_pol.gnu; gnuplot ../hist_pol_png.gnu"
    p = subprocess.Popen(cmdline, shell=True, cwd=wd)
    p.wait()
    
    cmdline = "mv hist_pol.png ../../plot/%s/%s_hist_pol.png" % (group,jobname)
    p = subprocess.Popen(cmdline, shell=True, cwd=wd)
    p.wait()
    
    cmdline = "mv hist_pol_1.png ../../plot/%s/%s_hist_pol_1.png" % (group,jobname)
    p = subprocess.Popen(cmdline, shell=True, cwd=wd)
    p.wait()

    cmdline = "mv hist_pol_2.png ../../plot/%s/%s_hist_pol_2.png" % (group,jobname)
    p = subprocess.Popen(cmdline, shell=True, cwd=wd)
    p.wait()
    
    cmdline = "mv hist_pol_3.png ../../plot/%s/%s_hist_pol_3.png" % (group,jobname)
    p = subprocess.Popen(cmdline, shell=True, cwd=wd)
    p.wait()
