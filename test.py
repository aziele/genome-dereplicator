#!/usr/bin/env python3

import pathlib
import unittest
import shutil
import tempfile

import gender as gd


class TestDistance(unittest.TestCase):

    def setUp(self):
        self.root_dir = pathlib.Path(__file__).parent
        self.in_dir = self.root_dir / 'test' / 'input'
        self.out_dir = pathlib.Path(tempfile.mkdtemp())

    def test_get_genome_paths(self):
        genome_paths = gd.get_genome_paths(self.in_dir)
        self.assertEqual(len(genome_paths), 10)

    def test_seqkit_stats(self):
        genome_paths = gd.get_genome_paths(self.in_dir)
        in_file = self.out_dir / 'seqkit.input.txt'
        out_file = self.out_dir / 'seqkit.stats.txt'
        process = gd.seqkit_stats(genome_paths, in_file, out_file, 1)
        self.assertEqual(process.returncode, 0)

    def test_get_n50s(self):
        genome_paths = gd.get_genome_paths(self.in_dir)
        in_file = self.out_dir / 'seqkit.input.txt'
        out_file = self.out_dir / 'seqkit.stats.txt'
        gd.seqkit_stats(genome_paths, in_file, out_file, 1)
        n50s = gd.get_n50s(out_file)
        self.assertEqual(n50s['1.fa'], 100000)
        self.assertEqual(n50s['2.fa'], 58150)
        self.assertEqual(n50s['4.fa'], 50000)

    def test_get_nonrepresentative_genomes(self):
        clusters = {'1': {'2', '3'}, '5': {'6'}, '8': set()}
        genomes = gd.get_nonrepresentative_genomes(clusters)
        self.assertEqual(genomes, {'2', '3', '6'})

    def test_update_clusters_1(self):
        prev = {'1': {'3'}, '6': set()}
        curr = {'4': set(), '7': {'6'}}
        update = gd.update_clusters(curr, prev)
        self.assertEqual(update, {'4': set(), '7': {'6'}, '1': {'3'}})

    def test_update_clusters_2(self):
        prev = {'2': {'3'}, '5': {'4'}}
        curr = {'1': {'2'}, '5': set()}
        update = gd.update_clusters(curr, prev)
        self.assertEqual(update, {'1': {'2', '3'}, '5': {'4'}})

    def test_update_clusters_3(self):
        prev = {}
        curr = {'1': {'2'}, '5': set()}
        update = gd.update_clusters(curr, prev)
        self.assertEqual(update, {'1': {'2'}, '5': set()})

    def test_kmerdb(self):
        genome_paths = gd.get_genome_paths(self.in_dir).values()
        db_file = self.out_dir / 'kmer.db'
        process = gd.kmerdb_build(genome_paths, db_file, 21, 1, self.out_dir, 1)
        self.assertEqual(process.returncode, 0)
        self.assertTrue(db_file.exists())
        out_file = self.out_dir / 'output.csv'
        process = gd.kmerdb_all2all(db_file, out_file, 1)
        self.assertEqual(process.returncode, 0)
        self.assertTrue(out_file.exists())
        process = gd.kmerdb_distance(out_file, 'ani', 0.95, 1)
        self.assertEqual(process.returncode, 0)

    def test_parse_distances(self):
        string = '''kmer-length: 21 fraction: 1 ,1.fa,2.fa,3.fa,4.fa,5.fa,6.fa,
        1.fa,
        2.fa,1:0.999998,
        3.fa,1:0.999986,2:0.999984,
        4.fa,
        5.fa,4:0.990897,
        6.fa,
        '''
        filename = self.out_dir / 'temp.csv'
        with open(filename, 'w') as oh:
            oh.write(string)
        distances = gd.parse_distances(filename)
        self.assertEqual(len(distances), 4)
        self.assertEqual(distances[0], ('2.fa', '1.fa', 0.999998))
        self.assertEqual(distances[1], ('3.fa', '1.fa', 0.999986))
        self.assertEqual(distances[2], ('3.fa', '2.fa', 0.999984))
        self.assertEqual(distances[3], ('5.fa', '4.fa', 0.990897))

    def test_create_graph(self):
        distances = [
            ('2.fa', '1.fa', 0.999998),
            ('3.fa', '1.fa', 0.999986),
            ('3.fa', '2.fa', 0.999984),
            ('5.fa', '4.fa', 0.990897)
        ]
        graph = {
            '2.fa': {'1.fa', '3.fa'},
            '1.fa': {'3.fa', '2.fa'},
            '3.fa': {'1.fa', '2.fa'},
            '5.fa': {'4.fa'},
            '4.fa': {'5.fa'}
        }
        test_graph = gd.create_graph(distances)
        self.assertEqual(graph, test_graph)

    def test_cluster_genomes(self):
        genomes = ['1.fa', '2.fa', '3.fa', '4.fa', '5.fa', '6.fa']
        distances = [
            ('2.fa', '1.fa', 0.999998),
            ('3.fa', '1.fa', 0.999986),
            ('3.fa', '2.fa', 0.999984),
            ('5.fa', '4.fa', 0.990897)
        ]
        graph = gd.create_graph(distances)
        clusters = [{'1.fa', '2.fa', '3.fa'}, {'4.fa', '5.fa'}, {'6.fa'}]
        self.assertEqual(clusters, gd.cluster_genomes(genomes, graph))

    def test_set_representative_genome_per_cluster(self):
        clusters = [{'2.fa', '1.fa', '3.fa'}, {'4.fa', '5.fa'}, {'6.fa'}]
        n50s = {'1.fa': 1000, '2.fa': 700, '3.fa': 800, '4.fa': 2000,
                '5.fa': 3000, '6.fa': 100}
        clusters = gd.set_representative_genome_per_cluster(clusters, n50s)
        self.assertEqual(clusters['1.fa'], {'3.fa', '2.fa'})
        self.assertEqual(clusters['5.fa'], {'4.fa'})
        self.assertEqual(clusters['6.fa'], set())

    def test_save_clusters(self):
        clusters = {
            'A.fa': {'C.fa', 'D.fa', 'B.fa'},
            'E.fa': set(),
            'F.fa': {'G.fa'}}
        n50s = {
            'A.fa': 4000, 'B.fa': 3000, 'C.fa': 2000, 'D.fa': 1000, 
            'E.fa': 9000, 'F.fa': 5000, 'G.fa': 4000}
        out_file = self.out_dir / 'clusters.txt'
        out_file = gd.save_clusters(clusters, out_file, n50s)
        with open(out_file) as fh:
            content = fh.read()
        self.assertEqual(content, 'A.fa*,B.fa,C.fa,D.fa\nE.fa*\nF.fa*,G.fa\n')

    def test_main_threshold_1(self):
        """ANI > 90%"""
        gd.main([
            '--ani',
            '0.9',
            '--disable-log-file',
            '--disable-log-stdout',
            f'{self.in_dir}',
            f'{self.out_dir}'
        ])
        repr_genomes = self.out_dir.joinpath('fasta').iterdir()
        repr_genomes = sorted([f.name for f in repr_genomes])
        repr_genomes_test = ['1.fa', '10.fa', '5.fa', '7.fa', '8.fa']
        self.assertEqual(repr_genomes, repr_genomes_test)

    def test_main_threshold_2(self):
        """ANI > 95%"""
        gd.main([
            '--ani',
            '0.95',
            '--disable-log-file',
            '--disable-log-stdout',
            f'{self.in_dir}',
            f'{self.out_dir}'
        ])
        repr_genomes = self.out_dir.joinpath('fasta').iterdir()
        repr_genomes = sorted([f.name for f in repr_genomes])
        repr_genomes_test = ['1.fa', '10.fa', '5.fa', '7.fa', '8.fa']
        self.assertEqual(repr_genomes, repr_genomes_test)

    def test_main_threshold_3(self):
        """ANI > 98%"""
        gd.main([
            '--ani',
            '0.98',
            '--disable-log-file',
            '--disable-log-stdout',
            f'{self.in_dir}',
            f'{self.out_dir}'
        ])
        repr_genomes = self.out_dir.joinpath('fasta').iterdir()
        repr_genomes = sorted([f.name for f in repr_genomes])
        repr_genomes_test = ['1.fa', '10.fa', '5.fa', '7.fa', '8.fa', '9.fa']
        self.assertEqual(repr_genomes, repr_genomes_test)

    def test_main_threshold_4(self):
        """ANI > 99%"""
        gd.main([
            '--ani',
            '0.99',
            '--disable-log-file',
            '--disable-log-stdout',
            f'{self.in_dir}',
            f'{self.out_dir}'
        ])
        repr_genomes = self.out_dir.joinpath('fasta').iterdir()
        repr_genomes = sorted([f.name for f in repr_genomes])
        test = ['1.fa', '10.fa', '5.fa', '6.fa', '7.fa', '8.fa', '9.fa']
        self.assertEqual(repr_genomes, test)

    def test_main_batch_1(self):
        gd.main([
            '--ani',
            '0.95',
            '--batch_size',
            '4',
            '--disable-log-file',
            '--disable-log-stdout',
            f'{self.in_dir}',
            f'{self.out_dir}'
        ])
        repr_genomes = self.out_dir.joinpath('fasta').iterdir()
        repr_genomes = sorted([f.name for f in repr_genomes])
        repr_genomes_test = ['1.fa', '10.fa', '5.fa', '7.fa', '8.fa']
        self.assertEqual(repr_genomes, repr_genomes_test)

    def test_main_batch_2(self):
        gd.main([
            '--ani',
            '0.95',
            '--batch_size',
            '6',
            '--disable-log-file',
            '--disable-log-stdout',
            f'{self.in_dir}',
            f'{self.out_dir}'
        ])
        repr_genomes = self.out_dir.joinpath('fasta').iterdir()
        repr_genomes = sorted([f.name for f in repr_genomes])
        repr_genomes_test = ['1.fa', '10.fa', '5.fa', '7.fa', '8.fa']
        self.assertEqual(repr_genomes, repr_genomes_test)

    def tearDown(self):
        shutil.rmtree(self.out_dir)


unittest.main()