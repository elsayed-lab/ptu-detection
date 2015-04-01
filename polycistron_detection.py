#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Trypanosome polycistron detection
Keith Hughitt <khughitt@umd.edu>
2013/06/01

The goal of this script is to assign each generate a classification scheme
which assigns each TriTryp gene to a category corresponding to its membership
in a polycistronic transcription unit.

For each chromosome, genes are clustered based on strand. Genes are then
assigned to tentative polycistronic units, or marked with 0 if it falls in
a gap region.
"""
import os
import pandas
import numpy as np
from scipy.ndimage.filters import median_filter

def main():
    """Main"""
    input_data = [
        {'species': 'T. brucei', 'file': 'TbruceiTreu927Gene_4.2_genes.csv'},
        {'species': 'T. cruzi',  'file': 'TcruziEsmeraldo_4.2_genes.csv'},
        {'species': 'L. major',  'file': 'LmajorFriedlinGene_4.2_genes.csv'}
    ]

    # specify window sizes to use
    window_sizes = [3, 5, 9, 15, 49, 99]

    # for each species...
    for input_ in input_data:
        filepath = os.path.join('../../data', input_['file'])

        genes = pandas.read_csv(filepath, sep='\t', skiprows=4)

        # create an empty data frame to append to
        cols = list(genes.columns) + ['med_%d' % i for i in window_sizes]
        output = pandas.DataFrame(columns=cols)

        # iterate through chromosomes and segment genes
        for i in genes.chromosome.unique():
            ch = genes[genes.chromosome == i]

            # segment genes in chromosome by strand
            ch['group'] = segment_data(ch.strand)

            # create a numeric version of the strand column
            strand = [1 if x=='+' else 0 for x in ch.strand]

            # use median filter with several window sizes in order to "smooth" the small
            # segments and create several alternate groupings with less fine scale structure
            for i in window_sizes:
                series = pandas.Series(median_filter(strand, i))
                ch['med_%d' % i] = np.asarray(segment_data(series))

            # break up groups that span gaps
            gap_size = 15000
            gap_indices = ch[ch.start.diff() > gap_size].index

            for idx in gap_indices:
                for i in window_sizes:
                    vals = (ch['med_%d' % i].loc[idx:] + 1).values
                    ch['med_%d' % i].loc[idx:] = vals
                ch.group.loc[idx:] = (ch.group.loc[idx:] + 1).values

            # append results to output data frame
            output = output.append(ch)

        # save result
        output_file = input_['file'].replace('_genes', '_segmented')
        output.to_csv(output_file, sep='\t')

def segment_data(x):
    """
    Returns an array representing the segments of similar values found in
    an input array.

    References
    ----------
    http://stackoverflow.com/questions/14358567/finding-consecutive-segments-in-a-pandas-data-frame
    """
    return (x.shift(1) != x).astype(int).cumsum()

if __name__ == "__main__":
    main()
    print("Done!")

