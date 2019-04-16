from mbf_genomics.genes.anno_tag_counts import (
    _FastTagCounter,
    IntervalStrategyExonSmart,
    IntervalStrategyExon,
    IntervalStrategyGene,
)

def get_interval_trees(interval_strategy, genome, chr):
    import bx.intervals

    by_chr = interval_strategy._get_interval_tuples_by_chr(genome)
    tree_forward = bx.intervals.IntervalTree()
    tree_reverse = bx.intervals.IntervalTree()
    gene_to_no = {}
    ii = 0
    for tup in by_chr[chr]:  # stable_id, strand, [starts], [stops]
        length = 0
        for start, stop in zip(tup[2], tup[3]):
            if tup[1] == 1:
                tree_forward.insert_interval(bx.intervals.Interval(start, stop, ii))
            else:
                tree_reverse.insert_interval(bx.intervals.Interval(start, stop, ii))
            length += stop - start
        gene_stable_id = tup[0]
        gene_to_no[gene_stable_id] = ii
        ii += 1
    return tree_forward, tree_reverse, gene_to_no

class _CounterStrategyBase:
    cores_needed = 1

    def count_reads(self, interval_strategy, genome, bamfile, reverse=False):
        lookup = {}
        for chr in genome.get_chromosome_lengths():
            lookup.update(
                self.count_gene_reads_on_chromosome(
                    interval_strategy, genome, bamfile, chr, False
                )
            )
        self.sanity_check(genome, lookup, bamfile, interval_strategy)
        return lookup

    def extract_lookup(self, data):
        """Adapter for count strategies that have different outputs
        (e.g. one-hashmap-unstranded or two-hashmaps-one-forward-one-reversed)
        """
        return data


class CounterStrategyStranded(_CounterStrategyBase):
    """This counter fetches() all reads on one chromosome at once, then matches them to the respective intervals
    defined by self.strategy.get_interval_trees"""

    def __init__(self):
        self.disable_sanity_check = False

    def count_gene_reads_on_chromosome(
        self, interval_strategy, genome, samfile, chr, reverse=False
    ):
        """Return a dict of gene_stable_id -> spliced exon matching read counts"""
        # basic algorithm: for each aligned region, check whether it overlaps a (merged) exon
        # if so, consider this read a hit for the gene with that exon
        # reads may hit multiple genes (think overlapping genes, where each could have generated the read)
        # multi aligned reads may count for multiple genes
        # but not for the same gene multiple times
        tree_forward, tree_reverse, gene_to_no = get_interval_trees(interval_strategy,
            genome, chr
        )

        counts = {}
        # one multi aligned read should only count once for a gene
        # further alignments into the same gene don't count
        # further alignments into other genes do though.
        # because of that' we need to keep the read names in a set.
        # uargh.
        for read in samfile.fetch(chr, 0, genome.get_chromosome_lengths()[chr]):
            seen = set()
            if not reverse:
                if read.is_reverse:
                    tree = tree_reverse
                else:
                    tree = tree_forward
            else:
                if read.is_reverse:
                    tree = tree_forward
                else:
                    tree = tree_reverse
            for sub_start, sub_stop in read.get_blocks():
                for x in tree.find(sub_start, sub_stop):
                    seen.add(
                        x.value
                    )  # one read matches one gene once. It might match multiple genes though.
            for ii in seen:
                if ii not in counts:
                    counts[ii] = set()
                counts[ii].add(read.qname)
                # counts[ii] += 1
        real_counts = {}
        no_to_gene = dict((v, k) for (k, v) in gene_to_no.items())
        for ii in counts:
            gene_stable_id = no_to_gene[ii]
            real_counts[gene_stable_id] = len(counts[ii])
        # if not real_counts and gene_to_no and len(chr) == 1:
        # raise ValueError("Nothing counted %s" % chr)
        return real_counts

    def sanity_check(self, genome, lookup, bam_file, interval_strategy):
        # the sanity check allows the exon tag count annotator to detect if you've reversed your reads
        if self.disable_sanity_check:
            return
        longest_chr = list(
            sorted([(v, k) for (k, v) in genome.get_chromosome_lengths().items()])
        )[-1][1]
        reverse_count = self.count_gene_reads_on_chromosome(
            interval_strategy, genome, bam_file, longest_chr, reverse=True
        )
        error_count = 0
        for gene_stable_id in reverse_count:
            if reverse_count[gene_stable_id] > 100 and reverse_count[
                gene_stable_id
            ] > 1.1 * (lookup[gene_stable_id] if gene_stable_id in lookup else 0):
                error_count += 1

        if error_count > 0.1 * len(reverse_count):
            # import cPickle
            # with open('debug.pickle','wb') as op:
            # cPickle.dump(lookup, op)
            # cPickle.dump(reverse_count, op)
            raise ValueError(
                "Found at least %.2f%% of genes on longest chromosome to have a reverse read count (%s) was above 110%% of the exon read count. This indicates that this lane should have been reversed before alignment. Set reverse_reads=True on your Lane object"
                % (100.0 * error_count / len(reverse_count), self.__class__.__name__)
            )


class CounterStrategyUnstranded(_CounterStrategyBase):
    """This counter fetches() all reads on one chromosome at once, 
    then matches them to the respective intervals
    defined by self.get_interval_trees"""

    def count_gene_reads_on_chromosome(
        self, interval_strategy, genome, samfile, chr, reverse=False
    ):
        """Return a dict of gene_stable_id -> spliced exon matching read counts"""
        # basic algorithm: for each aligned region, check whether it overlaps a (merged) exon
        # if so, consider this read a hit for the gene with that exon
        # reads may hit multiple genes (think overlapping genes, where each could have generated the read)
        # multi aligned reads may count for multiple genes
        # but not for the same gene multiple times
        tree_forward, tree_reverse, gene_to_no = get_interval_trees(interval_strategy,
            genome, chr
        )

        counts = {}
        # one multi aligned read should only count once for a gene
        # further alignments into the same gene don't count
        # further alignments into other genes do though.
        # because of that' we need to keep the read names in a set.
        # uargh.
        for read in samfile.fetch(chr, 0, genome.get_chromosome_lengths()[chr]):
            seen = set()
            for sub_start, sub_stop in read.get_blocks():
                for t in tree_forward, tree_reverse:
                    for x in t.find(sub_start, sub_stop):
                        seen.add(x.value)
            for ii in seen:
                if ii not in counts:
                    counts[ii] = set()
                counts[ii].add(read.qname)
                # counts[ii] += 1
        real_counts = {}
        no_to_gene = dict((v, k) for (k, v) in gene_to_no.items())
        for ii in counts:
            gene_stable_id = no_to_gene[ii]
            real_counts[gene_stable_id] = len(counts[ii])
        # if not real_counts and gene_to_no and len(chr) == 1:
        # print(counts)
        # print(real_counts)
        # print(not real_counts)
        #            raise ValueError("Nothing counted %s" % chr)
        return real_counts

    def sanity_check(self, genome, lookup, bam_file, interval_strategy):
        pass  # no op


class ExonSmartStrandedPython(_FastTagCounter):
    def __init__(self, aligned_lane):
        _FastTagCounter.__init__(
            self,
            aligned_lane,
            CounterStrategyStranded(),
            IntervalStrategyExonSmart(),
            "Exon, protein coding, stranded smart tag count %s",
            "Tag count inside exons of protein coding transcripts (all if no protein coding transcripts) exons, correct strand only",
        )


class ExonSmartUnstrandedPython(_FastTagCounter):
    def __init__(self, aligned_lane):
        _FastTagCounter.__init__(
            self,
            aligned_lane,
            CounterStrategyUnstranded(),
            IntervalStrategyExonSmart(),
            "Exon, protein coding, unstranded smart tag count %s",
            "Tag count inside exons of protein coding transcripts (all if no protein coding transcripts)  both strands",
        )


class ExonStrandedPython(_FastTagCounter):
    def __init__(self, aligned_lane):
        _FastTagCounter.__init__(
            self,
            aligned_lane,
            CounterStrategyStranded(),
            IntervalStrategyExon(),
            "Exon, protein coding, stranded tag count %s",
            "Tag count inside exons of protein coding transcripts (all if no protein coding transcripts) exons, correct strand only",
        )


class ExonUnstrandedPython(_FastTagCounter):
    def __init__(self, aligned_lane):
        _FastTagCounter.__init__(
            self,
            aligned_lane,
            CounterStrategyUnstranded(),
            IntervalStrategyExon(),
            "Exon, protein coding, unstranded tag count %s",
            "Tag count inside exons of protein coding transcripts (all if no protein coding transcripts)  both strands",
        )


class GeneStrandedPython(_FastTagCounter):
    def __init__(self, aligned_lane):
        _FastTagCounter.__init__(
            self,
            aligned_lane,
            CounterStrategyStranded(),
            IntervalStrategyGene(),
            "Gene, stranded tag count %s",
            "Tag count inside gene body (tss..tes), correct strand only",
        )


class GeneUnstrandedPython(_FastTagCounter):
    def __init__(self, aligned_lane):
        _FastTagCounter.__init__(
            self,
            aligned_lane,
            CounterStrategyUnstranded(),
            IntervalStrategyGene(),
            "Gene unstranded tag count %s",
            "Tag count inside gene body (tss..tes), both strands",
        )
