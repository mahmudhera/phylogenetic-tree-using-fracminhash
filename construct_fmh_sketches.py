import screed
import subprocess
import mmh3

class FracMinHash:
    '''
    FracMinHash class.
    '''
    def __init__(self, scale_factor, max_hash_value, initial_set=None):
        '''
        Create an FMH with given scale_facor in (0,1), largest hash value H.
        If initial_set is provided, that needs to be constructed as a set(),
        and must only contain integers in [0,H)
        Returns: None
        '''
        if initial_set is None:
            self.hash_set = set()
        else:
            self.hash_set = set(initial_set)
        self.H = max_hash_value
        self.scale_factor = scale_factor

    def add_value(self, hash_value):
        '''
        Add a hash value to the sketch.
        Returns: None
        '''
        if hash_value <= self.H * self.scale_factor:
            self.hash_set.add(hash_value)

    def add_values(self, hash_values):
        '''
        Add multiple hash values to the sketch.
        Returns: None
        '''
        for hash_value in hash_values:
            self.add_value(hash_value)

    def remove(self, hash_value):
        '''
        Remove a hash value from sketch.
        Returns: None
        '''
        self.hash_set -= hash_value

    def get_containment(self, smh):
        '''
        Obtain containment of provided sketch within this sketch.
        Returns: float -- containment value
        '''
        return 1.0 * len(self.hash_set.intersection(smh.hash_set)) / len(self.hash_set)

    def get_scaled_containment(self, smh, num_all_elements):
        '''
        Obtain scaled containment of provided sketch within this sketch.
        Accounts for the bias factor. Must provide L from the paper.
        Returns: float -- containment value
        '''
        bf = 1.0 - (1.0 - self.scale_factor) ** int(num_all_elements)
        return 1.0 * len(self.hash_set.intersection(smh.hash_set)) / ( len(self.hash_set) * bf )

    def get_sketch_size(self):
        '''
        Returns: int -- size of the sketch
        '''
        return len( self.hash_set )

def get_hash_from_kmer(kmer, seed=0):
    '''
    Using murmurhash, get hash value of a kmer and return that.
    H = largest hash value = 2^64
    Returns: int - a hash value
    Raises: ValueError if kmer has invalid characters
    '''
    kmer = kmer.upper()
    if not set(kmer).issubset('ATGC'):
        raise ValueError
    hash_value = mmh3.hash64(kmer, seed=seed)[0]
    if hash_value < 0:
        hash_value += 2**64
    return hash_value

def create_frac_minhash(kmers, seed, scale_facor):
    '''
    Given a list of kmers, generate the frac minhash sketch of those kmers.
    Returns: the sketch, an object of FracMinHash
    '''
    H = 2**64
    smh1 = FracMinHash(scale_facor, H)
    for kmer in kmers:
        h = get_hash_from_kmer(kmer, seed)
        smh1.add_value(h)
    return smh1

def add_kmers_in_scaled_minhash(kmers, smh, seed):
    '''
    Given a bunch of kmers and an fmh sketch, add these kmers to the sketch.
    Returns: the sketch, object of type FracMinHash
    '''
    H = 2**64
    for kmer in kmers:
        h = get_hash_from_kmer(kmer, seed)
        smh.add_value(h)
    return smh

def get_kmers_in_file_using_jellyfish(filename, k):
    # run jellyfish
    cmd = "jellyfish count -m " + str(k) + " -s 2G -t 32 -o tmp -C -L 1 -U 99999999 " + filename
    args = cmd.split(' ')
    subprocess.call(args)
    cmd = "jellyfish dump -L 1 -U 99999999 -o tmp-dump -c tmp"
    args = cmd.split(' ')
    subprocess.call(args)
    # open output file
    df = pd.read_csv('tmp-dump', delimiter=' ', header=None)
    list_kmers = df.iloc[:,0].tolist()
    return list_kmers

def build_kmers(sequence, ksize):
    '''
    Given a sequence and ksize, construct list of all kmers.
    Capital letters are used
    Skip all kmers with N or R.
    Returns: list - all kmers
    '''
    sequence = sequence.upper()
    kmers = []
    n_kmers = len(sequence) - ksize + 1
    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        if 'N' in kmer or 'R' in kmer:
            continue
        kmers.append(kmer)
    return kmers

def get_kmers_in_file_using_screed(filename, k=21):
    '''
    Given a filename and value of k, read all the kmers in the file and return as a list
    Returns: list - all kmers in file. Upper letters. N and R characters are skipped.
    '''
    all_kmers = []
    for record in screed.open(filename):
        sequence = record.sequence
        kmers = build_kmers(sequence, ksize)
        all_kmers += kmers
    return all_kmers

if __name__ == "__main__":
    seq1 = 'acgtgcgcgtgatgc'
    seq2 = 'acgttcgcgtgatgc'
    kmers1 = build_kmers(seq1, 2)
    kmers2 = build_kmers(seq2, 2)
    fmh1 = create_frac_minhash(kmers1, 1, 1.0)
    print(fmh1.hash_set)
    fmh2 = create_frac_minhash(kmers2, 1, 1.0)
    print(fmh2.hash_set)
    print( fmh1.get_containment(fmh2) )
    print( fmh2.get_containment(fmh1) )
