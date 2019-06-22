#!/usr/bin/python

import argparse
import re

def HistoIndex():
    parser = argparse.ArgumentParser()
    parser.add_argument("expr", nargs=1)
    args = parser.parse_args()
    expr = re.findall("[A-Za-z_0-9\[\]()\-<]+", args.expr[0])

    name = []
    total = []
    switch = False
    for item in expr:
        if item == "<":
            switch = True
            continue
        if not switch:
            result = re.findall("^[A-Za-z_(]+.*", item)
            if result:
                name += result
        else:
            result = re.findall("[0-9]+", item)
            if result:
                total += result

    length = len(name)
    if length != len(total):
        print("Length of name and total does not match")
        return

    ih = "int ih = "
    prod = ""
    for i in xrange(length):
        ih += prod + name[i]
        prod += total[i] + "*"
        if i < length-1:
            ih += " + "
        else:
            ih += " < " + prod[0:-1]
    print(ih)

if __name__ == "__main__":
    HistoIndex()
