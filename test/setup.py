from subprocess import call
from sys import argv

def findindex(a):
	return a[0]

call(['mpic++', '-g', 'testSetupCommunicators.cpp', '-o', 'split'])
with open('tmp.txt', 'w') as f:
	call(['mpirun', '-np', argv[1], './split', argv[2], argv[3], argv[4]], stdout=f)

with open('tmp.txt', 'r') as f:  # move cursor to beginning of file
	lines = f.readlines()
	lines = [line.split() for line in lines]
	lines = sorted(lines,key=findindex)
	for line in lines:
		print "worldRank:  " + line[0] + "  headRank:   " + line[1] + "  boundRank:" + line[2] + "  isFollower: " + line[3] + "  isLeader:   " + line[4]
call(['rm', 'tmp.txt'])
