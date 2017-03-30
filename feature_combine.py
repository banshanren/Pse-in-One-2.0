import argparse

parser = argparse.ArgumentParser(description = 'combine features from different files')
parser.add_argument('-f', help='gether features from one file or more', nargs='+')
args = parser.parse_args()
for i in args.f:
	print(i)