#! python

import os, sys, math, unittest
# from [library] import [function]
import dfstree

class test001(unittest.TestCase):
	def test_run(self):
		"""DFS initial test"""
		nodes = 5
		connections = [(1,0), (2,3), (3,4)]
		dfs = dfstree.DepthFirstSearchPartitioner()
		dfs.initialize(nodes)
		dfs.loadConnections(connections)
		sets = dfs.findConnectedSets()
		sets.sort(key=lambda x: len(x))
		self.assertTrue(len(sets[0]) == 2)
		self.assertTrue(len(sets[1]) == 3)
		self.assertTrue(sum(sets[0]) == 1)

	def test_run(self):
		"""DFS secondary test"""
		nodes = 10
		connections = [(5,6), (2,3), (3,4)]
		dfs = dfstree.DepthFirstSearchPartitioner()
		dfs.initialize(nodes)
		dfs.loadConnections(connections)
		sets = dfs.findConnectedSets()
		#print(sets)
		sets.sort(key=lambda x: len(x))
		self.assertTrue(len(sets) == 7)

	def test_nx(self):
		nodes = [x for x in range(100)]
		edges = [(1,2), (3,4), (1,5), (6,1)]
		subg = dfstree.findConnectedSubgraphs(nodes, edges)
		for s in subg:
			self.assertTrue(set(s) == set([3,4]) or set(s) == set([1,2,5,6]))
			#print(subg)
		#print([x for xs in subg for x in xs])

if __name__=="__main__":
	unittest.main(verbosity=2)