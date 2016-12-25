#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Trypanosome polycistron detection

Author:       Keith Hughitt <khughitt@umd.edu>
Created:      2013/06/01
Last updated: 2015/03/31

The goal of this script is to generate a classification scheme which assigns
each TriTryp gene to a category corresponding to its membership in a
polycistronic transcription unit.

For each chromosome, a 1d vector with length = # genes on the chromosome 
is created indicating the strands of genes as they appear on the chromosome,
moving from location 0 to the end of the chromosome. 

Next, smoothing is applied to add some robustness to small inter-PTU strand
differences (i.e. one to several genes on the opposing strand within a larger
cluster of same-strand genes).

Finally, the smoothed strand vectors are segmented into individual PTUs.

Based on some experimentation, a window-size of about 5 seems to product PTUs
which overlap well the actual underlying PTU's, as observed in IGV.

The input for this script are GFF annotations downloaded from TriTrypDB. Note
that if you are using an earlier version of TriTrypDB, the files may need to be
edited to remove any FASTA sequences appended below the GFF contents.

Example usage:

    python ptu_detection.py TriTrypDB-29_LmajorFriedlin_genes.gff

"""
import os
import sys
import pandas
import numpy as np
import warnings
from scipy.ndimage.filters import median_filter

def main():
    """Main"""
    # specify window sizes to use (unit = gene-strand)
    window_sizes = [1, 3, 5, 9, 15]

    # maximum gap size to allow within a single PTU
    gap_size = 10000

    # first argument should be filepath
    input_file = sys.argv[1]

    if not os.path.isfile(input_file):
        print("Invalid filepath specified")
        sys.exit()

    gff_fields = ['seqname', 'source', 'feature', 'start', 'end', 'score',
                  'strand', 'frame', 'attribute']
    gff = pandas.read_csv(input_file, sep='\t', comment='#', names=gff_fields)

    # tritrypdb organism abbreviation (e.g. LmjF)
    organism_id = gff.seqname[0].split('.').pop(0)

    # filter out non-genes
    #genes = gff[gff.attribute.str.startswith('ID=%s' % organism_id)]
    genes = gff[gff.feature == 'gene']

    # create an empty data frame to append to
    cols = ['gene_id'] + ['med_%d' % i for i in window_sizes]
    output = pandas.DataFrame(columns=cols)

    # create a list to keep track of BED file output rows
    bed_rows = []

    # iterate through chromosomes and segment genes
    for chr_name in genes.seqname.unique():
        chr_genes = genes[genes.seqname == chr_name]

        # resort by location (snoRNAs, etc. sometimes follow genes in GFF)
        chr_genes = chr_genes.sort('start').reset_index()

        # parse gene IDs
        with warnings.catch_warnings():
            gene_ids = []
            for attr_parts in chr_genes.attribute.str.split(';'):
                gene_ids.append(attr_parts[0][3:])
            
            # add gene_ids
            warnings.simplefilter("ignore")
            chr_genes['gene_id'] = gene_ids

            # segment genes in chromosome by strand
            warnings.simplefilter("ignore")
            chr_genes['group'] = segment_data(chr_genes.strand)

        # create a numeric version of the strand column
        strand = [1 if x=='+' else 0 for x in chr_genes.strand]

        # use median filter with several window sizes in order to "smooth" the small
        # segments and create several alternate groupings with less fine scale structure
        for win_size in window_sizes:
            series = pandas.Series(median_filter(strand, win_size))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                chr_genes['med_%d' % win_size] = np.asarray(segment_data(series))

        # break up groups that span gaps

        # Note: below code only compares CDS start positions so CDS will be
        # included in gap size calculation.
        #gap_indices = chr_genes[chr_genes.start.diff() > gap_size].index

        # gap_indices indicates index of gene immediately following a
        # sufficiently large gap
        intergenic_distances = chr_genes[1:].start.values - chr_genes[:-1].end.values
        is_gap = np.concatenate(([False], intergenic_distances > gap_size))
        gap_indices = chr_genes.index[is_gap]

        for idx in gap_indices:
            for i in window_sizes:
                vals = (chr_genes['med_%d' % i].loc[idx:] + 1).values
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    chr_genes['med_%d' % i].loc[idx:] = vals
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                chr_genes.group.loc[idx:] = (chr_genes.group.loc[idx:] + 1).values

        # add bed entries and rename PTUs
        for win_size in ["med_%d" % x for x in window_sizes]:
            # Create PTU name mapping
            ptu_mapping = {}
            ptu_counter = 1

            for ptu in chr_genes[win_size].unique():
                ptu_mapping[ptu] = "%s_PTU%02d" % (chr_name, ptu_counter)
                ptu_counter = ptu_counter + 1

            # Rename PTUs, e.g. "1.0" -> "LmjF.01_PTU01"
            #chr_genes[win_size] = chr_genes[win_size].apply(lambda x: "%s_PTU%02d" % (chr_name, x))
            chr_genes[win_size] = chr_genes[win_size].map(ptu_mapping)

            # Iterate over each PTU on the chromosome
            for ptu in chr_genes[win_size].unique():
                # Get genes for the PTU
                ptu_genes = chr_genes[chr_genes[win_size] == ptu]

                # Construct BED fields
                start_idx = str(int(ptu_genes.head(1).start))
                stop_idx = str(int(ptu_genes.tail(1).end))
                strand = ptu_genes.reset_index().head(1).strand[0]

                # Add bed entry
                bed_rows.append([
                    chr_name, start_idx, stop_idx, ptu, '0', strand
                ])

        # append results to output data frame
        output = output.append(chr_genes[cols])

    # save result
    output_file = os.path.basename(input_file).replace('.gff', '_PTUs.csv')

    # CSV
    output.to_csv(output_file, index=False)

    # Bed file
    fp = open(output_file.replace('.csv', '.bed'), 'w')

    fp.write('track name="Polycistronic Transcriptional Units" description=' +
             '"Polycistronic Transcriptional Units" ' +
             'visibility=2 colorByStrand="31,168,137 168,30,62"\n')
    for row in bed_rows:
        fp.write("\t".join(row) + "\n")

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

