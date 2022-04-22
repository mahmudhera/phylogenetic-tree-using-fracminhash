import screed
import subprocess
import mmh3


class FracMinHash:
    '''
    FracMinHash class.
    '''
    def __init__(self, scale_factor, max_hash_value, initial_set=None):
        '''
        Create an FMH with given scale_factor in (0,1), largest hash value H.
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

def create_frac_minhash(kmers, seed, scale_factor):
    '''
    Given a list of kmers, generate the frac minhash sketch of those kmers.
    Returns: the sketch, an object of FracMinHash
    '''
    H = 2**64
    smh1 = FracMinHash(scale_factor, H)
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

def get_kmers_in_file_using_screed(filename, ksize=21):
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

def read_genome_list(genome_list_filename):
    f = open(genome_list_filename, 'r')
    lines = f.readlines()
    num_genomes = int(lines[0].strip())
    genome_list = []
    for i in range(num_genomes):
        genome_name = lines[1 + i*2].strip()
        genome_path = lines[2 + i*2].strip()
        genome_list.append( (genome_name, genome_path) )
    f.close()
    return genome_list

def write_fmh_sketch(fmh, filename):
    '''
    Arguemnts: Object of type FracMinHash, and a filename
    Writes the hash values to that file
    Return None
    '''
    f = open(filename, 'w')
    f.write(str(fmh.scale_factor) + '\n')
    f.write(str(fmh.H) + '\n')
    for hash_value in fmh.hash_set:
        f.write(str(hash_value) + '\n')
    f.close()

def read_fmh_sketch(filename):
    '''
    Arguemnts: A filename
    Reads the hash values from that file
    Returs: an object of type FracMinHash
    '''
    hash_set = set()
    f = open(filename, 'r')
    lines = f.readlines()
    scale_factor = float( lines[0].strip() )
    max_hash_value = float( lines[1].strip() )
    for i in range(2, len(lines)):
        line = lines[i]
        hash_set.add( int( line.strip() ) )
    return FracMinHash(scale_factor, max_hash_value, hash_set)

def test_all():
    genome_name = 'test-genome/test_genome.fna'
    kmers = get_kmers_in_file_using_screed(genome_name, 21)
    print(len(kmers))
    fmh = create_frac_minhash(kmers, 0, 0.1)
    print(fmh.hash_set)
    write_fmh_sketch(fmh, 'test-genome/fmh-sketch1')
    fmh2 = read_fmh_sketch('test-genome/fmh-sketch1')
    print(fmh2.hash_set)
    print(fmh.hash_set.issubset( fmh2.hash_set ))
    print(fmh2.hash_set.issubset( fmh.hash_set ))

if __name__ == "__main__":
    genome_list_filename = 'genome-list-primates'
    sketch_directory = 'fmh_sketches'
    ksizes = [21, 31, 51]
    scale_factors = [0.000001, 0.00001]
    seeds = range(2)

    genome_list = read_genome_list(genome_list_filename)

    for (gname, gpath) in genome_list:
        f = open(gpath, 'r')
        print("File valid!")
        f.close()

    for (gname, gpath) in genome_list:
        print('Genome: ' + gname)
        count = 0
        for k in ksizes:
            for scale_factor in scale_factors:
                for seed in seeds:
                    scale = int(1.0/scale_factor)
                    sketch_filename = sketch_directory + '/fmh_sketch_k_' + str(k) + '_scale_f_' + str(scale_factor) + '_seed_' + str(seed) + '_genome_' + gname
                    cmd = 'sourmash compute -q -k '+ str(k) + ' --scaled ' + str(scale) + ' -o ' + sketch_filename + ' --seed ' + str(seed) + ' ' + gpath
                    print(cmd)
                    subprocess.call(cmd.split(' '))
        print("Done")
