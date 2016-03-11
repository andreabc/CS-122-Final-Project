from collections import defaultdict, Counter
from helpers.helpers import *
import cPickle as pickle
from os.path import join, exists, splitext
import time
import sys
import re

READ_LENGTH = 50

def hash_end(end, genome_ht):
    """
    Uses hashing to identify the set of locations spanned by
    a read.
    :param end: A single end of a read
    :param genome_ht: A hash of the genome with uniform key length
    :return:
    """
    key_length = len(genome_ht.keys()[0])
    end_pieces = [end[i * key_length: (i + 1) * key_length]
                  for i in range(len(end) / key_length)]

    hashed_read_locations = [genome_ht[read_piece]
                             for read_piece in end_pieces]
    start_positions = [[x - i * key_length for x in hashed_read_locations[i]]
                       for i in range(len(hashed_read_locations))]
    start_counter = Counter()

    for position_list in start_positions:
        start_counter.update(position_list)

    if not start_counter:
        return -1, 0
    else:
        best_alignment_location, best_alignment_count = \
            start_counter.most_common(1)[0]

    if best_alignment_count < 2:
        return -1, best_alignment_count
    else:
        return best_alignment_location, best_alignment_count


def hash_read(read, genome_ht):
    """
    Uses hashing to identify the set of locations spanned by
    a read.
    :param read: A single read
    :param genome_ht: A hash of the genome with uniform key length
    :return:
    """

    oriented_reads = [(read[0][::i], read[1][::j]) for i, j in ((1, -1), (-1, 1))]
    ## Either one end is forward and the other end is reversed, or vice versa.

    best_score = -1
    best_alignment_locations = (-1, -1)
    best_oriented_read = ('', '')
    for oriented_read in oriented_reads:
        hash_results = [hash_end(_, genome_ht) for _ in oriented_read]
        hash_locations = [_[0] for _ in hash_results]
        hash_score = sum([_[1] for _ in hash_results])
        if hash_score > best_score:
            best_alignment_locations = hash_locations
            best_oriented_read = oriented_read
    return best_oriented_read, best_alignment_locations


def make_genome_hash(reference, key_length):
    """
    :param reference: The reference as a string stored
    :param key_length: The length of keys to use.
    :return:
    """
    genome_hash = defaultdict(list)
    for i in range(len(reference) - key_length):
        ref_piece = reference[i: i + key_length]
        genome_hash[ref_piece].append(i)
    return genome_hash


def build_hash_and_pickle(ref_fn, key_length, force_rebuild=False):
    reference_hash_pkl_fn = '{}_hash.pkl'.format(splitext(ref_fn)[0])
    if exists(reference_hash_pkl_fn) and not force_rebuild:
        ref_genome_hash = pickle.load(open(reference_hash_pkl_fn, 'rb'))
        if len(ref_genome_hash.keys()[0]) == key_length:
            return ref_genome_hash
        else:
            pass
    else:
        pass
    reference = read_reference(ref_fn)
    ref_genome_hash = make_genome_hash(reference, key_length)
    pickle.dump(ref_genome_hash, open(reference_hash_pkl_fn, 'wb'))
    return ref_genome_hash


def hashing_algorithm(paired_end_reads, genome_ht, genome_length):
    alignments = []
    genome_aligned_reads = []
    count = 0
    start = time.clock()
    coverage = [0]*genome_length
    #print "coverage length = {}".format(len(coverage))

    for read in paired_end_reads:
        alignment, genome_aligned_read = hash_read(read, genome_ht)
        #print "alignment: {}, genome_aligned_read: {}".format(alignment, genome_aligned_read)
        if genome_aligned_read[0] != -1: 
            #coverage[genome_aligned_read[0]:genome_aligned_read[0]+READ_LENGTH+1] = [i + 1 for i in coverage[genome_aligned_read[0]:genome_aligned_read[0]+READ_LENGTH]]
            #print "{} {}".format(genome_aligned_read[0], genome_aligned_read[0] + 11)
            for i in range(genome_aligned_read[0], genome_aligned_read[0] + READ_LENGTH):
                if i < genome_length:
                    coverage[i] += 1
                #print i
        if genome_aligned_read[1] != -1:
            #coverage[genome_aligned_read[1]:genome_aligned_read[1]+READ_LENGTH+1] = [i + 1 for i in coverage[genome_aligned_read[0]:genome_aligned_read[0]+READ_LENGTH]]
            for i in range(genome_aligned_read[1], genome_aligned_read[1] + READ_LENGTH):
                if i < genome_length:
                    coverage[i] += 1
                #print i
        alignments.append(alignment)
        genome_aligned_reads.append(genome_aligned_read)
        count += 1
        if count % 100 == 0:
            time_passed = (time.clock()-start)/60
            print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print 'Approximately {:.3} minutes remaining'.format(remaining_time)
    
    for i in range(0,len(coverage),10):
        print i, coverage[i:i+10]
   
    return alignments, genome_aligned_reads, coverage, 

#ef repeats(CNV):
#   #2 bases
#   i = 0
#   start = CNV[0:2]
#   num_repeats = 0
#   while CNV[i:i+2] == start and i <= len(CNV)-2:
#       #print CNV[i:i+2]
#       i += 2
#       num_repeats += 1
#       if num_repeats > 4:
#           #STR_dict[start].append(CNV) 
#           return num_repeats, 2

#   #3 bases
#   i = 0        
#   start = CNV[0:3]
#   num_repeats = 0
#   while CNV[i:i+3] == start and i < len(CNV)-1:
#       #print CNV[i:i+3]
#       i += 3
#       num_repeats += 1
#       if num_repeats > 4:
#           #STR_dict[start].append(CNV)
#           return num_repeats, 3
#   #4 bases
#   i = 0
#   start = CNV[0:4]
#   num_repeats = 0
#   while CNV[i:i+4] == start and i < len(CNV)-1:
#       #print CNV[i:i+4]
#       i += 4
#       num_repeats += 1
#       if len(CNV) - num_repeats == 4:
#           #STR_dict[start].append(CNV)
#           return 4
#   #5 bases
#   i = 0
#   num_repeats = 0
#   start = CNV[0:5]
#   while CNV[i:i+5] == start and i < len(CNV)-1:
#       #print CNV[i:i+5]
#       i += 5
#       num_repeats += 1
#       if len(CNV) / num_repeats == 5:
#           #STR_dict[start].append(CNV)
#           return 5
#   return 0    

def is_STR(CNV):
    for x in range(len(CNV)):
        #2 bases
        i = 0
        start = CNV[x:x+2]
        num_repeats = 0
        while CNV[i:i+2] == start and i <= len(CNV)-2:
            #print CNV[i:i+2]
            i += 2
            num_repeats += 1
            if num_repeats > 4:
                #STR_dict[start].append(CNV) 
                return True

        #3 bases
        i = 0        
        start = CNV[x:x+3]
        num_repeats = 0
        while CNV[i:i+3] == start and i < len(CNV)-1:
            #print CNV[i:i+3]
            i += 3
            num_repeats += 1
            if num_repeats > 4:
                #STR_dict[start].append(CNV)
                return True
        #4 bases
        i = 0
        start = CNV[x:x+4]
        num_repeats = 0
        while CNV[i:i+4] == start and i < len(CNV)-1:
            #print CNV[i:i+4]
            i += 4
            num_repeats += 1
            if num_repeats > 4:
                #STR_dict[start].append(CNV)
                return True
        #5 bases
        i = 0
        num_repeats = 0
        start = CNV[x:x+5]
        while CNV[i:i+5] == start and i < len(CNV)-1:
            #print CNV[i:i+5]
            i += 5
            num_repeats += 1
            if num_repeats > 4:
                #STR_dict[start].append(CNV)
                return True
        

    return False
#def extend_STR(STR, donor_string, start_indices):
#    possible_new = STR
#    for x in start_indices:
#
#
#    return possible_new

def extend_CNV(CNV, num_matches, donor_string, start_indices):
    if len(start_indices) == 0:
        return CNV
    start_index = start_indices[0]

    possible_new = CNV
    end_index = start_index + len(CNV) -1

    #print "CNV = {}, start_index = {}, end_index = {}".format(CNV_or_STR, start_index, end_index)
    
    #extend left
    first_try = True
    str_start_indices = []
    j=1
    while (len(str_start_indices) >= num_matches or first_try) and start_index > 2:
        possible_new = donor_string[start_index-j] + CNV   #add on to the beginning
        p = re.compile(possible_new)
        it = re.finditer(p, donor_string)
        str_start_indices = [m.start() for m in it]
        #upstream = [donor_string[i-a:i] for i in str_start_indices]
        first_try = False
        j+=1
    if j == 2:  #if exited after first try
        possible_new = CNV

    #now extend right
    first_try = True
    str_start_indices = []

    j=1
    while (len(str_start_indices) >= num_matches or first_try) and end_index < len(donor_string) - 2:
        possible_new = donor_string[end_index+j] + CNV   #add on to the beginning
        p = re.compile(possible_new)
        it = re.finditer(p, donor_string)
        str_start_indices = [m.start() for m in it]
        first_try = False
        j+=1

    if j == 2:  #if exited after first try
        possible_new = CNV
    print "before: {} after: {}".format(CNV, possible_new)
    return possible_new

def get_CNV(ref_genome, ref_coverage, donor_genome):
    median_coverage = sorted(ref_coverage)[len(ref_coverage)/2]

    donor_string = ""
    with open(donor_genome, 'r') as f:
        for x in f:
            donor_string += x
    
    STR_dict = defaultdict(list)
    cnv_dict = defaultdict(list)

    for i in range(len(ref_coverage)-1):
        cnv = ""
        
        if ref_coverage[i] >= median_coverage*1.9:  #find middle of CNV
            #print "i = {}".format(i)
            cnv += ref_genome[i]
            j=i+1
            
            #continue right while median_coverage >= 1.5*median_coverage
            while ref_coverage[j] >= median_coverage*1.4 and j < len(ref_genome)-2:
                cnv+=ref_genome[j]
                j+=1
            
            #continue left of position
            j=i-1
            start = i
            while ref_coverage[j] >= median_coverage*1.4:
                cnv = ref_genome[j] + cnv
                j-=1
                start -= 1

            if len(cnv) > 20:
                cnv_dict[cnv].append(start)

        #if (ref_coverage[i+1]) - ref_coverage[i] > 3 or ref_coverage[i] > 1.7*median_coverage:           #find regions where there is a jump in coverage
        #    #print "position: {} - {}, {}".format(i+1,ref_coverage[i],ref_coverage[i+1])
        #    count = 0
        #    j = i+1
        #    while (ref_coverage[j] - ref_coverage[j+1] <= 3):   #find end of coverage increase
        #        if count > 800 or j == len(ref_coverage)-2:
        #            cnv = ""
        #            break
        #        cnv += ref_genome[j]
        #        count += 1
        #        j += 1
        #    if cnv != "" and len(cnv) > 20:
        #        cnv_dict[cnv].append(i+1)   #append start position to CNV dict
    
    a = 25  #upstream length
    
    for key in cnv_dict.keys():     #for every CNV
        #determine if extra coverage region is STR
        if is_STR(key):
            cnv_dict.pop(key)

            p = re.compile(key)
            it = re.finditer(p, donor_string)
            str_start_indices = [m.start() for m in it]

            #new_cnv = extend(key, len(str_start_indices), donor_string, str_start_indices[0])   #extend
            #p = re.compile(new_cnv)
            #it = re.finditer(p, donor_string)
            #str_start_indices = [m.start() for m in it]

            upstream = [donor_string[i-a:i] for i in str_start_indices]        

            for item in upstream:       #find upstream region in reference
                u = re.compile(item)    
                u_it = re.finditer(u, ref_genome)
                ref_upstream = [m.start() for m in u_it]
                ref_starts = [x + a for x in ref_upstream]  #add number to find STR ref position
                #print ref_starts

            if "ref_starts" in locals():
                for i in ref_starts:
                    STR_dict[key].append(i)
                STR_dict[key] = list(set(STR_dict[key]))


            continue
        
        #if CNV
        else:
            p = re.compile(key)          
            it = re.finditer(p,donor_string)
            cnv_start_indices = [m.start() for m in it]               #find all start positions of CNV in donor       
            
            new_cnv = extend_CNV(key, len(cnv_start_indices), donor_string, cnv_start_indices)   #extend
            p = re.compile(new_cnv)          
            it = re.finditer(p,donor_string)
            cnv_start_indices = [m.start() for m in it] 


            #upstream = [donor_string[i-a:i] for i in cnv_start_indices]    #get donor sequence 20bp ahead        
            #
            #for item in upstream:       #find upstream region in reference
            #    u = re.compile(item)    
            #    u_it = re.finditer(u, ref_genome)
            #    ref_upstream = [m.start() for m in u_it]
            #    ref_starts = [x + a for x in ref_upstream]  #add number to find CNV ref position
            #    #print ref_starts
            #first = cnv_dict[key][0]
            #cnv_dict.pop(key)
            #cnv_dict[key].append(first)
            #if "ref_starts" in locals():
            #    #print ref_starts
            #    for i in ref_starts:
            #        if i in cnv_dict[key]:
            #            continue
            #        cnv_dict[key].append(i)
            #        #q = re.compile(upstream)
            for i in cnv_start_indices:
                cnv_dict[key].append(i)
            cnv_dict[key] = list(set(cnv_dict[key]))
            
            if len(cnv_dict[key]) <= 1: #if CNV only appears in reference
                cnv_dict.pop(key)
            if len(key) < 40 and len(cnv_dict[key]) < 3:
                cnv_dict.pop(key)
            #if len(key) < 100 and len(cnv_dict[key]) < 2:
            #    cnv_dict.pop(key)

            #print "cnv_start_indices = {}".format(cnv_start_indices)
            #print cnv_dict[key]
            #if len(cnv_start_indices) <= 1: #if only found in ref
            #    cnv_dict.pop(key)
            #    continue
            #cnv_dict.pop(key)
            

    #print start_positions
    #print list_cnvs
    
    #for key in cnv_dict.keys():
    #    print key, cnv_dict[key]

    return cnv_dict, STR_dict

