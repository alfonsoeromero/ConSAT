import os
import unittest

from gfam.tasks.find_domain_arch.clustering_file import ClusteringFile


class TestClusteringFile(unittest.TestCase):
    def setUp(self):
        current_dir = os.path.dirname(__file__)
        data_dir = os.path.join(current_dir, os.pardir,
                                os.pardir, os.pardir, os.pardir, "data")
        self.cluster_file = os.path.join(data_dir, "cca.txt")

    def test_clustering_file_no_table_produces_expected_number(self):
        # arrange
        min_num_sequences: int = 3
        self._sut = ClusteringFile(min_num_sequences)
        unique_proteins = lambda x: {y.split(':')[0] for y in x.split()}
        expected_num_clusters = len([x for x in open(self.cluster_file)
                                     if len(unique_proteins(x)) >= min_num_sequences])

        # act
        data_table = self._sut.process_clustering_file(self.cluster_file)

        # assert
        self.assertEqual(expected_num_clusters, len(data_table))
