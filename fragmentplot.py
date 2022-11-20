#!/usr/bin/env python3
"""
# By Jana Alghoraibi
# This script parse alignment output in sam tab-delimited file format and generates a fragment recruitment-
plot using matplotlib graphing module. The script accept 3 arguments including: sam input file, reference genome name,
and sample id.
# Example on how to run the script:
python3 fragmentplot.py  SRX001355bowtie27.md.sam  'Parabacteroides distasonis' 'Sample: SRX001355'

"""

from matplotlib import pyplot as plt

import sys
import re

# Takes input sam file with aligned reads
input = sys.argv[1]
reference = sys.argv[2]
sampleid = sys.argv[3]

# Set lists for %identity and bp positions
y_match = []
x_match = []

genome = []

# Open file, read each line, extracts the alignment positions, the CIGAR and MD:Z columns to calculate % identity
with open(input) as input_file:
    # Parse each line in the file after skip headers
    for lines in input_file:

        if lines.startswith("@SQ"):
            l = re.search('(LN:)(\d+)', lines)
            genome_length = int(l.group(2))

            genome.append(genome_length)
            genome.sort()
            genome_size = genome[-1]

        if lines.startswith("@"):
            continue

        else:

            line = lines.rstrip("\n").split()

            # Extract mapping position
            position = int(line[3])

            x_match.append(position)

            # find CIGAR line and extract value of match/mismatch M

            CIGAR = line[5]
            m = re.search('([0-9]+)([M])', CIGAR)
            CIGAR_M = int(m.group(1))

            # Find MD:Z line and calculate the number of matches

            MD = line[13]
            n = re.search('[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*', MD)
            MDZ = n.group()
            matches = 0
            for entry in MDZ:
                if entry.isdigit():
                    matches += int(entry)

            identity = round((matches / CIGAR_M) * 100)

            # print(identity)

            y_match.append(identity)

# Create fragment recruitment plot using bp positions 'x-axis' and %identity 'y-axis'

plt.scatter(x_match, y_match)
plt.title(sampleid)
plt.suptitle(reference)
plt.xlabel('Genome Position (bp)')

plt.xlim(0, genome_size + 1)
plt.ylabel('% Identity')
plt.ylim([0, 101])

plt.show()
