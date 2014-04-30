from iva import kmers

class Contig:
    def __init__(self, fasta, verbose=0):
        self.fa = fasta
        self.verbose = verbose
        self.left_kmers = kmers.Kmers(left=True, verbose=verbose)
        self.right_kmers = kmers.Kmers(verbose=verbose)


    def __len__(self):
        return len(self.fa)


    def add_left_kmer(self, kmer):
        self.left_kmers.append(kmer)
        

    def add_right_kmer(self, kmer):
        self.right_kmers.append(kmer)


    def extend(self, min_cov, min_ratio, extend_bases):
        if self.verbose >= 2:
            print('    trying to extend left end ...')
        new_left_seq = self.left_kmers.extension(min_cov, min_ratio, extend_bases)
        if self.verbose >= 2:
            print('    trying to extend right end ...')
        new_right_seq = self.right_kmers.extension(min_cov, min_ratio, extend_bases)
        if self.verbose >= 1:
            print('    new left sequence:', new_left_seq)
            print('    new right sequence:', new_right_seq)
        self.fa.seq = new_left_seq + self.fa.seq + new_right_seq
        self.left_kmers = kmers.Kmers(left=True, verbose=self.verbose)
        self.right_kmers = kmers.Kmers(verbose=self.verbose)
        return len(new_left_seq), len(new_right_seq)
