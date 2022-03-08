#!/usr/bin/env python

import sys

from lop.file_io.ninfo import NinfoFile
from lop.elements.ninfo import NinfoSet


if len(sys.argv) != 3:
    print ("Usage: SCRIPT [input ninfo 1] [input ninfo 2]")
    sys.exit(2)


nf = NinfoFile(sys.argv[1])
nf.open_to_read()
ns1 = NinfoSet()
nf.read_all(ns1)
nf.close()

nf = NinfoFile(sys.argv[2])
nf.open_to_read()
ns2 = NinfoSet()
nf.read_all(ns2)
nf.close()

nhb1 = len(ns1.hbondDTs)
nhb2 = len(ns2.hbondDTs)

# Do not use index 0
flag1 = [False]*(nhb1+1)
flag2 = [False]*(nhb2+1)

def value_check(x1, x2):
    eps = 0.01
    if x1 - eps < x2 < x1 + eps:
        return True
    else:
        return False

for ihb1 in range(nhb1):

    hb1 = ns1.hbondDTs[ihb1]
    id1 = hb1.id
    pair1 = (hb1.imp1, hb1.imp2)

    for ihb2 in range(nhb2):

        hb2 = ns2.hbondDTs[ihb2]
        id2 = hb2.id

        if pair1 in ((hb2.imp1, hb2.imp2), (hb2.imp2, hb2.imp1)):

            if flag1[id1]:
                print('ihb1 is already found', id1)
                sys.exit(2)
            if flag2[id2]:
                print('ihb1 is already found', id2)
                sys.exit(2)

            # dist
            if not value_check(hb1.native, hb2.native):
                print('native differ', id1, id2, hb1.native, hb2.native)

            if not value_check(hb1.factor, hb2.factor):
                print('factor differ', id1, id2, hb1.factor, hb2.factor)

            if not value_check(hb1.coef, hb2.coef):
                print('coef differ', id1, id2, hb1.coef, hb2.coef)

            # angl
            # hb1.ang1 = hb2.ang1
            if (hb1.ang1_imp1, hb1.ang1_imp2, hb1.ang1_imp3) in ((hb2.ang1_imp1, hb2.ang1_imp2, hb2.ang1_imp3),
                                                                 (hb2.ang1_imp3, hb2.ang1_imp2, hb2.ang1_imp1)):

                if not value_check(hb1.ang1_native, hb2.ang1_native):
                    print('ang1_native differ', id1, id2, hb1.ang1_native, hb2.ang1_native)
                if not value_check(hb1.ang1_coef, hb2.ang1_coef):
                    print('ang1_coef differ', id1, id2, hb1.ang1_coef, hb2.ang1_coef)

                # hb1.ang2 must be the same as hb2.ang2
                if (hb1.ang2_imp1, hb1.ang2_imp2, hb1.ang2_imp3) not in ((hb2.ang2_imp1, hb2.ang2_imp2, hb2.ang2_imp3),
                                                                         (hb2.ang2_imp3, hb2.ang2_imp2, hb2.ang2_imp1)):
                    print('hb1.ang2 must equal to hb2.ang2')

                if not value_check(hb1.ang2_native, hb2.ang2_native):
                    print('ang2_native differ', id1, id2, hb1.ang2_native, hb2.ang2_native)
                if not value_check(hb1.ang2_coef, hb2.ang2_coef):
                    print('ang2_coef differ', id1, id2, hb1.ang2_coef, hb2.ang2_coef)


            # hb1.ang1 = hb2.ang2
            elif (hb1.ang1_imp1, hb1.ang1_imp2, hb1.ang1_imp3) in ((hb2.ang2_imp1, hb2.ang2_imp2, hb2.ang2_imp3),
                                                                   (hb2.ang2_imp3, hb2.ang2_imp2, hb2.ang2_imp1)):

                if not value_check(hb1.ang1_native, hb2.ang2_native):
                    print('ang1_native differ', id1, id2, hb1.ang1_native, hb2.ang2_native)
                if not value_check(hb1.ang1_coef, hb2.ang2_coef):
                    print('ang1_coef differ', id1, id2, hb1.ang1_coef, hb2.ang2_coef)

                # hb1.ang2 must be the same as hb2.ang1
                if (hb1.ang2_imp1, hb1.ang2_imp2, hb1.ang2_imp3) not in ((hb2.ang1_imp1, hb2.ang1_imp2, hb2.ang1_imp3),
                                                                         (hb2.ang1_imp3, hb2.ang1_imp2, hb2.ang1_imp1)):
                    print('hb1.ang2 must equal to hb2.ang1')

                if not value_check(hb1.ang2_native, hb2.ang1_native):
                    print('ang2_native differ', id1, id2, hb1.ang2_native, hb2.ang1_native)
                if not value_check(hb1.ang2_coef, hb2.ang1_coef):
                    print('ang2_coef differ', id1, id2, hb1.ang2_coef, hb2.ang1_coef)

            else:
                print("angl IDs do not match", id1, id2)

            # dihd
            # dih0
            if (hb1.dih0_imp1, hb1.dih0_imp2, hb1.dih0_imp3, hb1.dih0_imp4) not in ((hb2.dih0_imp1, hb2.dih0_imp2, hb2.dih0_imp3, hb2.dih0_imp4),
                                                                                    (hb2.dih0_imp4, hb2.dih0_imp3, hb2.dih0_imp2, hb2.dih0_imp1)):
                print('hb1.dih0 must equal to hb2.dih0')

            if not value_check(hb1.dih0_native, hb2.dih0_native):
                print('dih0_native differ', id1, id2, hb1.dih0_native, hb2.dih0_native)
            if not value_check(hb1.dih0_coef, hb2.dih0_coef):
                print('dih0_coef differ', id1, id2, hb1.dih0_coef, hb2.dih0_coef)

            # hb1.dih1 = hb2.dih1
            if (hb1.dih1_imp1, hb1.dih1_imp2, hb1.dih1_imp3, hb1.dih1_imp4) in ((hb2.dih1_imp1, hb2.dih1_imp2, hb2.dih1_imp3, hb2.dih1_imp4),
                                                                                (hb2.dih1_imp4, hb2.dih1_imp3, hb2.dih1_imp2, hb2.dih1_imp1)):

                if not value_check(hb1.dih1_native, hb2.dih1_native):
                    print('dih1_native differ', id1, id2, hb1.dih1_native, hb2.dih1_native)
                if not value_check(hb1.dih1_coef, hb2.dih1_coef):
                    print('dih1_coef differ', id1, id2, hb1.dih1_coef, hb2.dih1_coef)

                # hb1.dih2 must equal to hb2.dih2
                if (hb1.dih2_imp1, hb1.dih2_imp2, hb1.dih2_imp3, hb1.dih2_imp4) not in ((hb2.dih2_imp1, hb2.dih2_imp2, hb2.dih2_imp3, hb2.dih2_imp4),
                                                                                        (hb2.dih2_imp4, hb2.dih2_imp3, hb2.dih2_imp2, hb2.dih2_imp1)):
                    print('hb1.dih2 must equal to hb2.dih2')

                if not value_check(hb1.dih2_native, hb2.dih2_native):
                    print('dih2_native differ', id1, id2, hb1.dih2_native, hb2.dih2_native)
                if not value_check(hb1.dih2_coef, hb2.dih2_coef):
                    print('dih2_coef differ', id1, id2, hb1.dih2_coef, hb2.dih2_coef)

            # hb1.dih1 = hb2.dih2
            elif (hb1.dih1_imp1, hb1.dih1_imp2, hb1.dih1_imp3, hb1.dih1_imp4) in ((hb2.dih2_imp1, hb2.dih2_imp2, hb2.dih2_imp3, hb2.dih2_imp4),
                                                                                  (hb2.dih2_imp4, hb2.dih2_imp3, hb2.dih2_imp2, hb2.dih2_imp1)):

                if not value_check(hb1.dih1_native, hb2.dih2_native):
                    print('dih1_native differ', id1, id2, hb1.dih1_native, hb2.dih2_native)
                if not value_check(hb1.dih1_coef, hb2.dih2_coef):
                    print('dih1_coef differ', id1, id2, hb1.dih1_coef, hb2.dih2_coef)

                # hb1.dih2 must equal to hb2.dih1
                if (hb1.dih2_imp1, hb1.dih2_imp2, hb1.dih2_imp3, hb1.dih2_imp4) not in ((hb2.dih1_imp1, hb2.dih1_imp2, hb2.dih1_imp3, hb2.dih1_imp4),
                                                                                        (hb2.dih1_imp4, hb2.dih1_imp3, hb2.dih1_imp2, hb2.dih1_imp1)):
                    print('hb1.dih2 must equal to hb2.dih1')

                if not value_check(hb1.dih2_native, hb2.dih1_native):
                    print('dih2_native differ', id1, id2, hb1.dih2_native, hb2.dih1_native)
                if not value_check(hb1.dih2_coef, hb2.dih1_coef):
                    print('dih2_coef differ', id1, id2, hb1.dih2_coef, hb2.dih1_coef)

            else:
                print("dih IDs do not match", id1, id2)


            flag1[id1] = True
            flag2[id2] = True




for id1 in range(1, nhb1+1):
    if flag1[id1]:
        pass
    else:
        print('No match found in ninfo2, but exists in ninfo1:', id1)

for id2 in range(1, nhb2+1):
    if flag2[id2]:
        pass
    else:
        print('No match found in ninfo1, but exists in ninfo2:', id2)
