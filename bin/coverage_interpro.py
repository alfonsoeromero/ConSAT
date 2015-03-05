import sys
from collections import defaultdict
from collections import Counter

def fasta(f):
    sequences = []
    id = ""
    for line in f:
        if line.startswith(">"):
            new_id = line.strip().split("|")[1]
            if sequences:
                yield id, "".join(sequences)
            sequences = []
            id = new_id
        elif line.strip():
            sequences.append(line.strip())
    if sequences:
        yield id, "".join(sequences)

def merge_ranges(ranges):
    """
    Function taken from http://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap-in-python
    Merge overlapping and adjacent ranges and yield the merged ranges
    in order. The argument must be an iterable of pairs (start, stop).

    >>> list(merge_ranges([(5,7), (3,5), (-1,3)]))
    [(-1, 7)]
    >>> list(merge_ranges([(5,6), (3,4), (1,2)]))
    [(1, 2), (3, 4), (5, 6)]
    >>> list(merge_ranges([]))
    []
    """
    ranges = iter(sorted(ranges))
    current_start, current_stop = next(ranges)
    for start, stop in ranges:
        if start > current_stop:
            # Gap between segments: output current segment and start a new one.
            yield current_start, current_stop
            current_start, current_stop = start, stop
        else:
            # Segments adjacent or overlapping: merge.
            current_stop = max(current_stop, stop)
    yield current_start, current_stop

def count_residues(assignments):
    return sum([end - start + 1 for start, end in merge_ranges(assignments)])         

def extract_assignments(f):
    prev = ""
    prev_len = 0
    assignments = []
    for line in f:
        fields = line.strip().split("\t", 9)
        prot, source, start, end = fields[0], fields[3], int(fields[6]), int(fields[7])
        if prev and prot != prev:
            yield prev, assignments
            assignments = []
        assignments.append((source, start, end))
        prev = prot
    if assignments:
        yield prev, assignments

def count_architectures(output_file, arch_table_file, sequences, residues):
    proteins_per_arch = defaultdict(int)
    sequences_covered = 0
    residues_covered = 0

    sequences_covered_gfams = 0
    archs_with_gfams = set()
    archs_with_only_gfams = set()
    covered_residues_gfams = 0
    sequences_only_gfams = 0
    
    for line in open(arch_table_file):
        fields = line.strip().split("\t")
        len_covered, arch, arch_with_assignments = int(fields[2]), fields[3], fields[5]
        if arch == "NO_ASSIGNMENT":
            continue
        sequences_covered += 1
        if sequences_covered % 100000 == 0:
            print "Sequences covered by ConSAT processed: {}".format(sequences_covered)
        residues_covered += len_covered

        proteins_per_arch[arch] += 1
        if "CPPD" in arch:
            archs_with_gfams.add(arch)
            sequences_covered_gfams += 1
            only_gfam = True
            for domain in arch_with_assignments.replace('}', ' ').replace('{', ' ').replace(';', ' ').split():
                if domain.startswith("CPPD"):
                    start, end = domain.split("(")[1].replace(')', '').split('-')
                    covered_residues_gfams += (int(end) - int(start) + 1)
                else:
                    only_gfam = False
            if only_gfam:
                archs_with_only_gfams.add(arch)
                sequences_only_gfams += 1
    
    
    with open(output_file, "a") as out:
        out.write("\n\n")
        out.write("Number of sequences covered by ConSAT: {} ({}%)\n".format(sequences_covered, float(sequences_covered)*100.0/sequences))
        out.write("Number of residues covered by ConSAT: {} ({}%)\n".format(residues_covered, float(residues_covered)*100.0/residues))
        out.write("\n\n")
        out.write("Number of different architectures found: {}\n".format(len(proteins_per_arch)))
        out.write("\t of which {} ({}%) are participated by CPPDs domains\n".format(len(archs_with_gfams), float(len(archs_with_gfams))/len(proteins_per_arch)*100.0))
        out.write("\t where {} are architectures only formed by CPPDs\n".format(len(archs_with_only_gfams)))
        out.write("\n\n")
        out.write("CPPDs domains cover {} residues ({}%)\n".format(covered_residues_gfams, float(covered_residues_gfams)*100.0/residues))
        out.write("\t from {} sequences ({}%)\n".format(sequences_covered_gfams, float(sequences_covered_gfams)*100.0/sequences))
        out.write("\t of which {} are covered only with CPPDs\n".format(sequences_only_gfams))

    c = Counter(proteins_per_arch.values())
    with open("values", "w") as out:
        for value, cnt in c.items():
            out.write("{}\t{}\n".format(value, cnt))
    #plt.plot(*zip(*sorted(c.items()))
    #plt.savefig("./architectures_plot.png")
    #plt.close()

def count(ipr_file, output_file, sequences, residues):
    sources = set()
    sequences_with_assignments = 0
    residues_with_assignments = 0
    sequences_per_source = defaultdict(int)
    residues_per_source = defaultdict(int)
        
    f_interpro = open(ipr_file)
    for sequences_with_assignments, (protein, assignments) in enumerate(extract_assignments(f_interpro)):
        for source in {x[0] for x in assignments}:
            # one more sequence for this source
            sequences_per_source[source] += 1
            # we extract all the associated intervals
            intervals_source = [(start, end) for s, start, end in assignments if s == source]
            # we sum to the total residues per source
            residues_per_source[source] += count_residues(intervals_source)
        # we count all covered residues
        all_intervals = [(start, end) for s, start, end in assignments]
        residues_with_assignments += count_residues(all_intervals)
        if sequences_with_assignments % 100000 == 0:
            print "Number of sequences in IPR file processed: {}".format(sequences_with_assignments)

    with open(output_file, "w") as out:
        out.write("Total sequences: {}, residues: {}\n".format(sequences, residues))
        out.write("Total sequences with assignments: {} ({}%)\n".format(sequences_with_assignments, sequences_with_assignments/float(sequences)*100.0))
        out.write("Total residues covered: {} ({}%)\n".format(residues_with_assignments, residues_with_assignments/float(residues)*100.0))
        for source, num in sequences_per_source.items():
            out.write("Total sequences with assignments from source {}: {} ({}%)\n".format(source, num, float(num)/sequences*100.0))
            out.write("Total residues covered by source {}: {}\n".format(source, residues_per_source[source], float(residues_per_source[source])*100.0/float(residues)))

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print "USE: {} fasta_file_uniprot intepro_file architectures_table output_file".format(sys.argv[0])
        sys.exit(-1)
    fasta_file, interpro_file, arch_table, output_file = sys.argv[1:]

    # 1.- count sequences and residues
    fi = open(fasta_file)
    residues = 0
    for sequences, (id, seq) in enumerate(fasta(fi)):
        residues += len(seq)
        sequences += 1
        if sequences % 100000 == 0:
            print "Num processed sequences: {}".format(sequences)
    print "Residues: {}, sequences: {}".format(residues, sequences)

    # 2.- count coverage of IPR
    count(interpro_file, output_file, sequences, residues)
    for line in open(output_file):
        print line.strip()
    
    # 3.- loot at the architectures
    count_architectures(output_file, arch_table, sequences, residues)
    print
    for line in open(output_file):
        print line.strip()
