#!/usr/bin/env python3
"""Genome Dereplicator (gender) is a tool for dereplicating genomes."""

from __future__ import annotations
import argparse
import collections
import logging
import multiprocessing
import pathlib
import random
import shutil
import subprocess
import sys


APP_NAME = 'gender'
__version__ = '0.0.2'


def get_arguments(args):
    desc = f'{APP_NAME} v.{__version__}'
    parser = argparse.ArgumentParser(description=desc, add_help=False)

    p = parser.add_argument_group('Positional arguments (required)')
    p.add_argument('in_dir',
                   help='Input directory with genome assemblies')
    p.add_argument('out_dir', 
                   help='Output directory (will be created if not exists)')

    p = parser.add_argument_group('Optional arguments')
    p.add_argument('--ani', type=float, default=0.95,
                   help='ANI threshold (Average Nucleotide Identity). It is a '
                   'minimum ANI to cluster genomes [default: %(default)s]')
    p.add_argument('--k', type=int, default=21,
                   help='k-mer size [default: %(default)s]')
    p.add_argument('--k_fraction', type=float, default=1,
                   help='Fraction of all k-mers to include in dereplication '
                   '[default: %(default)s]')    
    p.add_argument('--measure', choices=['mash', 'min', 'jaccard'],
                   default='mash',
                   help='Measure to estimate ANI [default: %(default)s]')
    p.add_argument('--batch_size', type=int, default=50000,
                   help='Iterated dereplication on random batches of this many '
                   'genomes. Smaller number will reduce memory usage')         
    p.add_argument('--threads', type=int, 
                   default=min(multiprocessing.cpu_count(), 16),
                   help='Number of CPU threads [default: %(default)s]')

    p = parser.add_argument_group('Other arguments')
    p.add_argument('--disable-log-stdout', action='store_false',
                   help='Disable logging to stdout (quiet)')
    p.add_argument('--disable-log-file', action='store_false',
                   help='Disable logging to a file')
    p.add_argument('--debug', action='store_true',
                   help='Print debug messages to stdout')
    p.add_argument('-h', '--help', action='help',
                   default=argparse.SUPPRESS,
                   help='Print this help message and exit')
    p.add_argument('-v', '--version', action='version',
                   version=f'v{__version__}',
                   help="Print version information and exit")

    args = parser.parse_args(args)
    return args


def create_logger(
        name: str,
        file_name: Union[str, pathlib.Path, False],
        stdout: Optional[bool] = True,
        log_level: int = logging.INFO
        ) -> logging.Logger:
    """Returns a logger to log events.

    Args:
        name:
            Name of the logger.
        file_name:
            Name/Path of the logfile to be created. If the argument evaluates
            to False, logging to a file will be disabled.
        stdout:
            True/False enables/disables logging to stdout.
        log_level:
            The numeric level of the logging event (one of DEBUG, INFO etc.).

    """
    logger = logging.getLogger(name)
    logger.setLevel(log_level)

    # Set log format to handlers
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')

    if file_name:
        # Create file logger handler
        fh = logging.FileHandler(file_name, 'w')
        fh.setLevel(log_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    if stdout:
        # Create stream logger handler
        sh = logging.StreamHandler()
        sh.setLevel(log_level)
        sh.setFormatter(formatter)
        logger.addHandler(sh)

    return logger


def dereplicate(
        genome_list: list[str],
        genome_paths: dict[str, pathlib.Path],
        genome_n50s: dict[str, float],
        threshold: float,
        measure: Literal['mash', 'jaccard', 'min'],
        k_size: int, 
        k_fraction: float,
        temp_dir: pathlib.Path,
        threads: int
        ) -> dict[str, set[str]]:
    """Returns genome clusters with the representative genome in each cluster.

    Args:
        genome_list:
            A list of genome names to dereplicate (cluster).
        genome_paths:
            A dict mapping genome names to the corresponding file path.
        genome_n50s:
            A dict mapping genome names to the corresponding n50 scores.
        threshold:
            Minimum sequence identity (ANI) to cluster genomes (0-1].
        measure:
            Method to estimate ANI.
        k_size:
            K-mer size.
        k_fraction:
            Fraction of k-mers to include in dereplication (0, 1].
        temp_dir:
            Directory in which temporary files will be created.
        threads:
            Number of CPU threads.

    Returns:
        A dict mapping representative genomes to the corresponding cluster. Each
        cluster is represented as a set of genome names. For example:

        {'genome1': {'genome2', 'genome3'},
         'genome4': {'genome5'},
         'genome6': set()}
    """
    logger = logging.getLogger(APP_NAME)
    measure = 'ani' if measure == 'mash' else measure

    # Kmer-db build
    db_file = temp_dir.joinpath('kmer.db')
    genome_paths = [genome_paths[g] for g in genome_list]
    logger.debug('Running kmer-db build')
    process = kmerdb_build(
        genome_paths,
        db_file,
        k_size,
        k_fraction, 
        temp_dir,
        threads)
    if process.returncode:
        logger.error(process.stderr)
        sys.exit(1)

    # Kmer-db run
    out_file = temp_dir.joinpath('output.csv')
    logger.debug('Running kmer-db all2all')
    process = kmerdb_all2all(db_file, out_file, threads)
    if process.returncode:
        logger.error(process.stderr)
        sys.exit(1)

    # Kmer-db dist
    logger.debug('Running kmer-db distance')
    process = kmerdb_distance(out_file, measure, threshold, threads)
    out_file = pathlib.Path(f'{out_file}.{measure}')
    if process.returncode:
        logger.error(process.stderr)
        sys.exit(1)

    # Parse distances from the kmer-db output file
    logger.debug('Loading distances')
    pairwise_distances = parse_distances(out_file)
    # Create graph from distances
    logger.debug('Creating a graph from distances')
    graph = create_graph(pairwise_distances)
    # Cluster genomes
    logger.debug('Clustering genomes')
    clusters = cluster_genomes(genome_list, graph)
    # Set a representative genome in each cluster
    logger.debug('Setting a representative genome per cluster')
    clusters = set_representative_genome_per_cluster(clusters, genome_n50s)
    return clusters


def set_representative_genome_per_cluster(
        cluster_list: list[set[str]],
        genome_n50s: dict[str, float]
        ) -> dict[str, set[str]]:
    """Sets a representative genome in every cluster based on the largest N50.
    
    Args:
        cluster_list:
            A list of clusters. Each cluster is a set of genome names.
        genome_n50s:
            A dict mapping genome names to the corresponding n50 scores.

    Returns:
        A dict mapping representative genomes to the corresponding cluster. Each
        cluster is represented as a set of genome names. For example:

        {'genome1': {'genome2', 'genome3'},
         'genome4': {'genome5'},
         'genome6': set()}
    """
    clusters = {}
    for genomes in cluster_list:
        genomes = list(genomes)
        if len(genomes) > 1:
            genomes.sort(key=lambda x: genome_n50s[x], reverse=True)
        else:
            assert len(genomes) == 1
        clusters[genomes[0]] = set(genomes[1:])
    return clusters


def update_clusters(curr: dict[str, set[str]], prev: dict[str, set[str]]):
    """Updates current clusters based on clusters from previous dereplication.
    
    The function is used in dereplication running on random batches of genomes.
    In this case, the function updates clusters from a current batch with 
    clusters from the previous batch. Thus, the function mutates one of its 
    arguments (`curr`). In the case of dereplication occurring on all genomes, 
    the function has no side effect on any of its arguments.

    Args:
        curr:
            Current cluster.
        prev:
            Previous cluster.

    Returns:
        A dict mapping representative genomes to the corresponding cluster.
        Each cluster is represented as a set of genome names.
    """
    if not prev: return curr
    for repr_genome in curr:
        for name in list(curr[repr_genome]):
            if name in prev:
                curr[repr_genome].update(prev[name])
                del prev[name]
    for repr_genome in prev:
        if repr_genome not in curr:
            curr[repr_genome] = set()
        curr[repr_genome].update(prev[repr_genome])
    return curr


def get_nonrepresentative_genomes(clusters: dict[str, set[str]]) -> set(str):
    """Returns a set of all non-representative genomes present in the clusters.

    Note that non-representative genomes were already clustered (dereplicated)
    and thus they will be excluded from further dereplication iterations.
    """
    s = set()
    for genome_set in clusters.values():
        s.update(genome_set)
    return s


def save_clusters(
        clusters: dict[str, set[str]],
        out_path: pathlib.Path,
        genome_n50s: dict[str, float]
        ) -> pathlib.Path:
    """Saves clusters to a text file."""
    with open(out_path, 'w') as oh:
        for genome in clusters:
            oh.write(f'{genome}*')  # A representative genome in the cluster
            # Sort other genomes in the cluster by N50
            lst = sorted(clusters[genome], key=lambda x: genome_n50s[x], 
                reverse=True)
            if clusters[genome]:
                oh.write(f',{",".join(lst)}')
            oh.write('\n')
    return out_path


def kmerdb_build(
        genome_paths: list[str],
        db_file: pathlib.Path,
        k_size: int,
        k_fraction: float,
        temp_dir: pathlib.Path,
        threads: int
        ) -> subprocess.CompletedProcess:
    """Runs Kmer-db build."""
    in_path = temp_dir.joinpath('input.txt')
    with open(in_path, 'w') as oh:
        for file_path in genome_paths:
            oh.write(f'{file_path}\n')
    cmd = [
        'kmer-db',
        'build',
        '-f',
        f'{k_fraction}',
        '-k',
        f'{k_size}',
        '-t',
        f'{threads}',
        f'{in_path}',
        f'{db_file}'
    ]
    return subprocess.run(cmd, capture_output=True, text=True)


def kmerdb_all2all(
        db_file: pathlib.Path,
        out_file: pathlib.Path,
        threads: int
        ) -> subprocess.CompletedProcess:
    """Runs Kmer-db all2all."""
    cmd = [
        'kmer-db',
        'all2all',
        '-sparse',
        '-t',
        f'{threads}',
        db_file,
        out_file
    ] 
    return subprocess.run(cmd, capture_output=True, text=True)


def kmerdb_distance(
        in_file: pathlib.Path,
        measure: str,
        threshold: float,
        threads: int
        ) -> subprocess.CompletedProcess:
    """Runs Kmer-db distance."""
    cmd = [
        'kmer-db',
        'distance',
        '-sparse',
        '-above',
        f'{threshold}',
        measure,
        '-t',
        f'{threads}',
        f'{in_file}'
    ]
    return subprocess.run(cmd, capture_output=True, text=True)


def parse_distances(in_file: filepath.Path) -> list[tuple[str, str, float]]:
    """Returns a list of paiwise distances between genomes.

    Args:
        in_file:
            A path to the distance file from Kmer-db.
    """
    lst = []
    with open(in_file) as fh:
        ids = fh.readline().rstrip(',').split(',')[1:]
        idx2id = {i:id for i, id in enumerate(ids)}
        for line in fh:
            cols = line.strip().rstrip(',').split(',')
            id1 = cols[0].strip()
            for i, col in enumerate(cols[1:]):
                scol = col.split(':')
                idx = int(scol[0])
                val = float(scol[1])
                id2 = idx2id[idx-1]
                lst.append((id1, id2, val))
    return lst


def get_genome_paths(in_dir: pathlib.Path) -> dict[str, pathlib.Path]:
    """Finds genome FASTA files in a given directory.

    Args:
        in_dir:
            Input directory

    Returns:
        A dict mapping genome names to the corresponding file path. For example:

        {'1.fa': PosixPath('test/input/1.fa'),
         '2.fna': PosixPath('test/input/2.fa')}
    """
    formats = ['.fasta', '.fasta.gz', '.fa', '.fa.gz', '.fna', '.fna.gz']
    d = {}
    for file_path in sorted(in_dir.glob('**/*')):
        if file_path.is_file():
            for frmt in formats:
                if file_path.name.endswith(frmt):
                    d[file_path.name] = file_path
                    break
    return d


def seqkit_stats(
        genome_paths: dict[str, pathlib.Path],
        in_file: pathlib.Path,
        out_file: pathlib.Path,
        threads: int
        ) -> subprocess.CompletedProcess:
    """Runs seqkit stats on input genome files

    Args:
        genome_paths:
            A dict mapping genome names to the corresponding file path.
        in_file:
            Ppath to input file (in seqkit: --infile-list)
        out_file:
            Path to output file
        threads:
            Number of CPU threads.
    """
    with open(in_file, 'w') as oh:
        for file_name, file_path in genome_paths.items():
            oh.write(f'{file_path}\n')
    cmd = [
        'seqkit',
        'stats',
        '--all',
        '--infile-list',
        f'{in_file}',
        '--threads',
        f'{threads}',
        '--seq-type',
        'dna',
        '--tabular',
        '-o',
        f'{out_file}'
    ]
    process = subprocess.run(cmd, capture_output=True, text=True)
    return process


def get_n50s(in_file: pathlib.Path) -> dict[str, float]:
    """Returns N50s for input genomes based on seqkit output.

    Args:
        in_file:
            Path to seqkit output.
    
    Returns:
        A dict mapping genome names to the corresponding n50 scores.
        For example:

        {'1.fa': 1000000, '2.fna': 50000}
    """
    # Parse N50s from seqkit output
    d = {}
    fh = open(in_file)
    fh.readline()  # Skip header
    for line in fh:
        cols = line.split()
        path = pathlib.Path(cols[0])
        n50 = float(cols[12])
        d[path.name] = n50
    fh.close()
    return d


def create_graph(
        pairwise_distances: list[tuple[str, str, float]]
        ) -> dict[str, set[str]]:
    """Creates a graph of genomes which have identity above the threshold.

    An unidirected graph where nodes are genomes and edges connect genomes.

    Args:
        pairwise_distances:
            A list of paiwise distances between genomes.

    Returns:
        A dict mapping a given genome to a set of similar genomes (above the 
        threshold to that genome). For example:

        pairwise_distances = [
          ('2.fa', '1.fa', 0.999998),
          ('3.fa', '1.fa', 0.999986),
          ('3.fa', '2.fa', 0.999984),
          ('5.fa', '4.fa', 0.990897)
        ]
        graph = create_graph(pairwise_distances)
        print(graph)
        {'2.fa': {'1.fa', '3.fa'},
         '1.fa': {'3.fa', '2.fa'},
         '3.fa': {'1.fa', '2.fa'},
         '5.fa': {'4.fa'},
         '4.fa': {'5.fa'}} 
    """
    graph = collections.defaultdict(set)
    for g1, g2, distance in pairwise_distances:
        graph[g1].add(g2)
        graph[g2].add(g1)
    return graph


def cluster_genomes(
        genomes: list[str],
        graph: dict[str, set[str]]
        ) -> list[set[str]]:
    """Clusters genomes based on a given graph via single-linkage clustering.

    Args:
        genomes:
            A list of genomes.
        graph:
            A dict mapping a given genome to a set of similar genomes.

    Returns:
        A list of clusters. Each cluster is a set of genomes. For example:

        [{'1.fa', '2.fa', '3.fa'}, {'4.fa', '5.fa'}, {'6.fa'}]
    """
    visited = set()
    clusters = []
    for genome in genomes:
        if genome in visited:
            continue
        connected = dfs(graph, genome)
        clusters.append(connected)
        visited |= connected
    return clusters


def dfs(graph: dict[str, set[str]], start: str) -> set[str]:
    """Depth-Fist Search to visit all nodes."""
    visited, stack = set(), {start}
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.update(graph[vertex] - visited)
    return visited


def main(args=None):
    args = get_arguments(args)
    random.seed(0)
    
    # Set paths to working directories (input, output, and temp)
    in_dir = pathlib.Path(args.in_dir)
    out_dir = pathlib.Path(args.out_dir)
    out_subdir = out_dir / 'fasta'
    out_subdir.mkdir(parents=True, exist_ok=True)
    temp_dir = out_dir / 'temp'
    temp_dir.mkdir(parents=True, exist_ok=True)

    # Set log options
    log_file = out_dir / 'log.txt'
    logger = create_logger(
        APP_NAME,
        log_file if args.disable_log_file else False,
        args.disable_log_stdout,
        logging.DEBUG if args.debug else logging.INFO)

    # Log the user parameters
    params = [
        f'ANI: {args.ani}',
        f'k: {args.k}',
        f'k_fraction: {args.k_fraction}',
        f'batch_size: {args.batch_size}',
        f'threads: {args.threads}'
    ]
    logger.info(', '.join(params))

    # Get paths to input genome files
    ALL_GENOME_PATHS = get_genome_paths(in_dir)
    ALL_GENOME_COUNT = len(ALL_GENOME_PATHS)
    logger.info(f'Input genomes: {len(ALL_GENOME_PATHS):,}')

    # Calculate N50s for input genome files
    logger.info('Calculating N50s for input genomes')
    in_file = temp_dir / 'seqkit.input.txt'
    out_file = out_dir / 'seqkit.stats.txt'
    process = seqkit_stats(ALL_GENOME_PATHS, in_file, out_file, args.threads)
    if process.returncode:
        logger.error(process.stderr)
        sys.exit(1)
    logger.debug('Loading N50s')
    n50s = get_n50s(out_file)

    genome_list = list(ALL_GENOME_PATHS.keys())
    excluded_genomes = set()
    prev_clusters = {}
    # The following while loop is used only for dereplication in the batch mode.
    # Briefly, its goal is to cluster genomes within smaller randomly-selected
    # batches. This is repeated iteratively until an iteration fails to remove 
    # any genomes. When this happens, the script assumes that most replication 
    # has been cleared out and conducts a final dereplication round with all 
    # remaining genomes.
    while True:
        if len(genome_list) <= args.batch_size:
            break

        logger.info(f'Clustering {args.batch_size} genomes [random batch]')
        random.shuffle(genome_list)
        batch_genomes = genome_list[:args.batch_size]
        curr_clusters = dereplicate(
            batch_genomes,
            ALL_GENOME_PATHS,
            n50s,
            args.ani,
            args.measure,
            args.k,
            args.k_fraction,
            temp_dir,
            args.threads)
        
        if not (nonrep_genomes := get_nonrepresentative_genomes(curr_clusters)):
            logger.info(f'No genome clusters found.')
            break
        else:
            excluded_genomes |= nonrep_genomes

        genome_list = [x for x in genome_list if x not in excluded_genomes]
        logger.info(f'Remaining genomes: {len(genome_list):,}')
        # Update clusters
        curr_clusters = update_clusters(curr_clusters, prev_clusters)
        prev_clusters = curr_clusters

    # Final dereplication
    logger.info(f'Clustering {len(genome_list):,} genomes')
    curr_clusters = dereplicate(
        genome_list, 
        ALL_GENOME_PATHS,
        n50s,
        args.ani,
        args.measure,
        args.k,
        args.k_fraction,
        temp_dir, 
        args.threads)
    excluded_genomes |= get_nonrepresentative_genomes(curr_clusters)
    genome_list = [x for x in genome_list if x not in excluded_genomes]
    curr_clusters = update_clusters(curr_clusters, prev_clusters)
    logger.info(
        f'Clusters: {len(genome_list):,} ({ALL_GENOME_COUNT:,} genomes)')

    # Save clusters to file
    out_file = out_dir.joinpath('clusters.txt')
    logger.info(f'Saving clusters to {out_file}')
    save_clusters(curr_clusters, out_file, n50s)

    # Copy representative genomes
    logger.info(f'Copying representative genome per cluster to {out_subdir}')
    for genome in genome_list:
        shutil.copy(ALL_GENOME_PATHS[genome], out_subdir)

    # Remove temp_dir
    shutil.rmtree(temp_dir)

    if args.disable_log_file:
        logger.info(f'Log file: {log_file}')


if __name__ == '__main__':
    main()
