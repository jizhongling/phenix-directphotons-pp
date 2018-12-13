#!/usr/bin/python

import scipy.stats as stats

def RawALL():
    f = open("data/raw-asym.txt", "r")

    ipt = 0
    while True:
        lines = []
        count = 0

        while count < 8:
            line = f.readline().split(',')
            if len(line) <= 1:
                return
            line.pop()
            line = [float(item) for item in line]
            lines.append(line)
            count += 1

        print( ipt, stats.f_oneway(lines[0], lines[1], lines[2], lines[3],
              lines[4], lines[5], lines[6], lines[7]) )
        ipt += 1

if __name__ == "__main__":
    RawALL()
