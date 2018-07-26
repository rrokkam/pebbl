from subprocess import call
from sys import argv

def findindex(a):
	return a[0]

call(['mpic++', '-g', 'setup.cpp', '-o', 'split'])
with open('tmp.txt', 'w') as f:
	call(['mpirun', '-np', argv[1], './split', argv[2], argv[3], argv[4], argv[5]], stdout=f)

with open('tmp.txt', 'r') as f:  # move cursor to beginning of file
	lines = f.readlines()
	lines = [map(int, line.split()) for line in lines]
	lines = sorted(lines,key=findindex)
	for line in lines:
		print "worldRank: " + str(line[0]) + "  headRank: " + str(line[1])\
		+ "  boundRank: " + str(line[2]) + "  isFollower: " + str(line[3]) + \
		"  isLeader: " + str(line[4])
call(['rm', 'tmp.txt'])

