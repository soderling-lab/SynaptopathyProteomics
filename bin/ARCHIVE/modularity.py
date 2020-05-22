#!/usr/bin/env python3

# Python wrapper around the ada executable for modularity calculation.

## Usage:
#  $ ./modularity_calculation.exe network.net clusters.clu type

from argparse import ArgumentParser
ap = ArgumentParser(description=''' ''')
ap.add_argument("network", help = ''' ''')
ap.add_argument("clusters", help = ''' ''')
ap.add_argument("type", help = ''' ''')

args = vars(ap.parse_args())
data = args['data']

cmd = []
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
out = process.communicate()
