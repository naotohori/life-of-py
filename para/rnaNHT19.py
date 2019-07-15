class NHT19(object):
    # U_BL
    BL_PS = 23.0
    BL_SB = 10.0
    BL_SP = 64.0

    # U_BA
    BA_PSB =  5.0
    BA_PSP = 20.0
    BA_BSP =  5.0
    BA_SPS = 20.0

    # U_ST (consecutive)
    ST_DIST = 1.4
    ST_DIH = 4.0
    #*** st_u0 is defined in another place. See <<<< NHT19_stack_param

    # U_ST (non-consecutive)
    TST_DIST = 5.00
    TST_ANGL = 1.50
    TST_DIH = 0.15
    TST_U0 = -5.00

    # U_HB
    HB_DIST = 5.0
    HB_ANGL = 1.5
    HB_DIH_HBOND = 0.15
    HB_DIH_CHAIN = 0.15
    HB_U0 = -2.7
    #hb_cutoff_dist = 2.0


#<<<< NHT19_stack_param 
#**** 
#**** U0 = - h + kB ( T - Tm ) * s
#**** 
#**    h        s       Tm
#AA  3.92529  -0.319  298.913
#AC  3.89356  -0.319  298.913
#AG  4.68538   5.301  341.183
#AU  3.89356  -0.319  298.913
#CA  3.84872  -0.319  298.913
#CC  3.58574  -1.567  285.829
#CG  4.16163   0.774  315.519
#CU  3.55578  -1.567  285.829
#GA  4.64425   5.301  341.183
#GC  4.64896   4.370  343.196
#GG  5.12502   7.346  366.344
#GU  4.54742   2.924  338.164
#UA  3.84872  -0.319  298.913
#UC  3.54836   1.567  285.829
#UG  4.60496   2.924  338.164
#UU  2.92091  -3.563  251.610
#>>>>
#
#<<<< NHT19_exv_param
#**      R      epsilon
#P     1.8900   0.200000
#S     2.6100   0.200000
#A     2.5200   0.200000
#G     2.7000   0.200000
#C     2.4300   0.200000
#U     2.4300   0.200000
#Mg2   2.0000   0.894700
#Ca2   2.8000   1.000000
#CHX   5.0000   1.000000
#X1    2.1     0.2
#>>>>
