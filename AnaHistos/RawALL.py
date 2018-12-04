import scipy.stats as stats

f = open("data/raw-asym.txt", "r")
lines = []

while True:
  line = f.readline().split(',')
  if len(line) <= 1:
    break
  line.pop()
  line = [float(item) for item in line]
  lines.append(line)

for ipt in xrange(0,11):
  print( ipt, stats.f_oneway(lines[ipt*8+0], lines[ipt*8+1], lines[ipt*8+2], lines[ipt*8+3],
        lines[ipt*8+4], lines[ipt*8+5], lines[ipt*8+6], lines[ipt*8+7]) )
