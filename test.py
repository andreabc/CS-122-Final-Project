
def is_STR(CNV, STR_dict):
    #2 bases
    i = 0
    num_repeats = 0
    while CNV[i:i+2] == STR_dict[i] and i <= len(CNV)-2:
        start = CNV[0:2]
        i += 2
        num_repeats += 1
        if num_repeats > 3:
        	STR_dict[start].append(CNV) 
        	return True

    #3 bases
    i = 0
    num_repeats = 0
    while CNV[i:i+3] == STR_dict[i] and i < len(CNV)-1:
        start = CNV[0:3]
        i += 3
        num_repeats += 1
        if num_repeats > 3:
        	STR_dict[start].append(CNV)
        	return True
    #4 bases
    i = 0
    num_repeats = 0
    while CNV[i:i+4] == STR_dict[i] and i < len(CNV)-1:
        start = CNV[0:4]
        i += 4
        num_repeats += 1
        if num_repeats > 3:
        	STR_dict[start].append(CNV)
        	return True

    #5 bases
    i = 0
    num_repeats = 0
    while CNV[i:i+5] == STR_dict[i] and i < len(CNV)-1:
        start = CNV[0:5]
        i += 5
        num_repeats += 1
        if num_repeats > 3:
        	STR_dict[start].append(CNV)
        	return True

    return False

if __name__ == '__main__':
	STR_dict = defaultdict(list)
	STR_file = open("possible_SNPS.txt", 'r'):   
		for line in STR_file:
            STR_dict[line.strip()] = ""

    CNV1 = "ATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCA"
    CNV2 = "CAATCGTGTTCATTCCAGCCACAACTGACAACTCA"

    if is_STR(CNV1, STR_dict):
    	print "CNV1 is an STR"
    else:
    	print "CNV1 is not an STR"
    if is_STR(CNV2, STR_dict):
    	print "CNV1 is an STR"
    else: 
    	print "CNV2 is not an STR"