import sys
import numpy as np
from os.path import join
import time
from collections import defaultdict

def read_reads(read_fn):
    f = open(read_fn, 'r')
    first_line = True
    all_reads = []
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()     #strips whitespace
        paired_end_reads = line.split(',')  # The two paired ends are separated by a comma
        all_reads.append(paired_end_reads)
    return all_reads


def read_reference(ref_fn):     #concatenates reference
    f = open(ref_fn, 'r')
    first_line = True
    output_reference = ''
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()
        output_reference += line  # We append each line to the output reference string.
    return output_reference


def trivial_algorithm(paired_end_reads, ref):
    """

    This is a functional aligner, but it's a huge simplification that
    generate a LOT of potential bugs.  It's also very slow.

    Read the spec carefully; consider how the paired-end reads are
    generated, and ideally, write your own algorithm
    instead of trying to tweak this one (which isn't very good).

    :param paired_end_reads: Paired-end reads generated from read_reads
    :param ref: A reference genome generated from read_reference
    :return: 2 lists:
                1) a list of alignment locations for each read (all_alignment_locations).
                    The list gives the starting position of the minimum-mismatch alignment of both reads.
                2) a list of the paired-end reads set so that both reads are in their optimal orientation
                   with respect to the reference genome.
    """
    #hash all positions of reference (50bp)
#    count = 0
#    ref_dict = defaultdict(float)
#    for i in range(len(ref)):
#        ref_frag = ""
#        count_fifty = 0
#        if (count >= (len(ref)-50)):
#            break
#        while count_fifty <= 50:
#            ref_frag += ref[count + count_fifty]
#            count_fifty += 1
#        ref_dict[ref_frag] = count
#        count += 1
#
#    all_read_alignment_locations = []
#    output_read_pairs = []
#    count = 0
#    start = time.clock()
#    for read_pair in paired_end_reads:
#        count += 1
#        read_alignment_locations = []
#        output_read_pair = []
#        if count % 20 == 0:
#            time_passed = (time.clock()-start)/60
#            print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)
#            remaining_time = time_passed/count*(len(paired_end_reads)-count)
#            print 'Approximately {:.2} minutes remaining'.format(remaining_time
#        for read in read_pair:
#            if ref_dict[read]:
#                all_read_alignment_locations.append(ref_dict[read])
#
#            reversed_read = read[::-1]
#            if (ref_dict[reversed_read]):
#                all_read_alignment_locations.append(ref_dict[reversed_read])




    all_read_alignment_locations = []
    output_read_pairs = []
    count = 0
    start = time.clock()
    for read_pair in paired_end_reads:
        count += 1
        read_alignment_locations = []
        output_read_pair = []
        if count % 20 == 0:
            time_passed = (time.clock()-start)/60
            print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print 'Approximately {:.2} minutes remaining'.format(remaining_time)
        for read in read_pair:
            min_mismatches = 50
            min_mismatch_location = -1
            for i in range(len(ref) - 50):
                mismatches = [1 if read[j] != ref[i + j] else 0 for j in range(len(read))]
                n_mismatches = sum(mismatches)
                # The above line should be familiar to Python users, but bears  some explanation for
                # people who are getting started with it. The "mismatches = ..." line
                # is called a "list comprehension. Basically, this is a short way of writing the loop:
                #
                # n_mismatches = 0
                # for j in range(len(read)):
                # if read[j] != ref[i+j]:
                #         n_mismatches += 1
                #
                # The first line creates a list which has a 1 for every mismatch and a 0 for every match.
                # The second line sums the list created by the first line, which counts the number of mismatches.
                if n_mismatches < min_mismatches:
                    min_mismatches = n_mismatches
                    min_mismatch_location = i

            reversed_read = read[::-1]
            for i in range(len(ref) - 50):
                mismatches = [1 if reversed_read[j] != ref[i + j] else 0 for j in range(len(read))]
                n_mismatches = sum(mismatches)
                if n_mismatches < min_mismatches:
                    min_mismatches = n_mismatches
                    min_mismatch_location = i
                    read = reversed_read
            read_alignment_locations.append(min_mismatch_location)
            output_read_pair.append(read)
            # # Note that there are some huge potential problems here.

        all_read_alignment_locations.append(read_alignment_locations)
        output_read_pairs.append(output_read_pair)
    return all_read_alignment_locations, output_read_pairs


def pretty_print_aligned_reads_with_ref(genome_oriented_reads, read_alignments, ref):
    """

    :param genome_oriented_reads: oriented reads generated by trivial_algorithm
    :param read_alignments: alignments generated from trivial_algorithm
    :param ref: reference generated by read_ref
    :return: Returns nothing, but prints the reads aligned to the genome to
     show you what pileup actually *LOOKS* like. You should be able to call SNPs
     by eyeballing the output. However, there are some reads that will not align.
     In the future you'll want to re-check why these reads aren't aligning--the cause
     is usually a structural variation, like an insertion or deletion.
    """
    output_str = ''
    good_alignments = [120 < x[1] - x[0] < 180 for x in read_alignments]
    # There should be 50 + x (90 < x < 110) p between the reads, and we give a little
    # extra space in case there's been a deletion or insertion.  Depending on the type of
    # deletions/insertions

    best_reads = [genome_oriented_reads[i] for i in range(len(good_alignments))
                  if good_alignments[i]]
    # Remove the reads that do not have a good alignment, or a good reverse alignment.
    best_alignments = [read_alignments[i] for i in range(len(read_alignments))
                       if good_alignments[i]]
    # Take their corresponding alignments
    aligned_reads = [best_reads[i][0] + '.' * (best_alignments[i][1] - best_alignments[i][0] - 50)
                     + best_reads[i][1] for i in range(len(best_reads))]
    # This turns the reads into strings oriented towards the genome.
    # We get the first read, followed by the correct number of dots to join the first and second reads,
    # and then the second read.

    first_alignment = [x[0] for x in best_alignments]
    alignment_indices = np.argsort(first_alignment)
    sorted_reads = [aligned_reads[i] for i in alignment_indices]
    sorted_alignments = [best_alignments[i] for i in alignment_indices]

    # You don't need to worry too much about how the code block below works--its job is to make it so
    # that a read that starts printing in the third row will continue printing in the third row of the
    # next set of lines.
    active_reads = []
    line_length = 100
    output_str += '\n\n' + '-' * (line_length + 6) + '\n\n'
    for i in range(len(ref) / line_length):
        next_ref = ref[i * line_length: (i + 1) * line_length]
        new_read_indices = [j for j in range(len(sorted_reads))
                            if i * line_length <= sorted_alignments[j][0] < (i + 1) * line_length]
        space_amounts = [sorted_alignments[index][0] % line_length for index in new_read_indices]
        new_reads = [sorted_reads[index] for index in new_read_indices]
        new_reads_with_spaces = [' ' * space_amounts[j] + new_reads[j] for j in range(len(new_reads))]
        empty_active_read_indices = [index for index in range(len(active_reads)) if active_reads[index] == '']
        for j in range(min(len(new_reads_with_spaces), len(empty_active_read_indices))):
            active_reads[empty_active_read_indices[j]] = new_reads_with_spaces[j]

        if len(new_reads_with_spaces) > len(empty_active_read_indices):
            active_reads += new_reads_with_spaces[len(empty_active_read_indices):]
        printed_reads = ['Read: ' + read[:line_length] for read in active_reads]
        active_reads = [read[line_length:] for read in active_reads]
        while len(active_reads) > 0:
            last_thing = active_reads.pop()
            if last_thing != '':
                active_reads.append(last_thing)
                break
        output_lines = ['Ref:  ' + next_ref] + printed_reads
        output_str += 'Reference index: ' + str(i * line_length) + \
                      '\n' + '\n'.join(output_lines) + '\n\n' + '-' * (line_length + 6) + '\n\n'
    # print output_str
    return output_str


if __name__ == "__main__":
    folder = 'practice'
    f_base = '{}_W_chr_1'.format(folder)  #.format(folder) goes inside curly braces
    reads_fn = join(folder, 'reads_{}.txt'.format(f_base))
    start = time.clock()
    input_reads = read_reads(reads_fn)
    # This is for speed;
    # If you want to read everything
    # remove the [:300] part of the above line.

    reference_fn = join(folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
    # donor_fn = join(folder, 'donor_{}'.format(f_base))
    # donor = read_reference(donor_fn)
    trivial_algorithm(input_reads, reference)
    #alignments, reads = trivial_algorithm(input_reads, reference)
    #print alignments
    #print reads
    #output_str = pretty_print_aligned_reads_with_ref(reads, alignments, reference)
    #output_fn = join(folder, 'aligned_{}.txt'.format(f_base))
    ## print output_fn
    #with(open(output_fn, 'w')) as output_file:
    #    output_file.write(output_str)