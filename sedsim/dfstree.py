import sys
from scipy.special import factorial
import networkx as nx

class Node(object):
	def __init__(self, id):
		self.id = id
		self.connections = []
		self.visited = False

	def add(self, id):
		self.connections.append(id)

	def __str__(self):
		return "{:d}".format(self.id)

# Implements the algorithm from
# https://www.geeksforgeeks.org/connected-components-in-an-undirected-graph/
class DepthFirstSearchPartitioner(object):
	def __init__(self):
		self.nodes = {}
		self._counter = 0

	def initialize(self, num_nodes):
		for i in range(num_nodes):
			n = Node(i)
			self.nodes[i] = n

	def loadConnections(self, pair_list):
		ordered_pair_list = []
		for (i,j) in pair_list:
			self.nodes[i].add(j)
			self.nodes[j].add(i)

	def markRecursive(self, node, connections):
		node.visited = True
		self._counter += 1
		connections.append(node)
		for nid in node.connections:
			n = self.nodes[nid] 
			if not n.visited:
				#print("Recursively visiting node {:d} ({:d})".format(n.id, self._counter))
				self.markRecursive(n, connections)

	def findConnectedSets(self):
		# Required to prevent recursion errors
		#print("recur in")
		lim = sys.getrecursionlimit()
		if lim < len(self.nodes):
			print(lim, len(self.nodes))
			#sys.setrecursionlimit(len(self.nodes))
		#print(len(self.nodes))
		connected_ids = []
		for nid in self.nodes:
			node = self.nodes[nid]
			connections = []
			if not node.visited:
				self.markRecursive(node, connections)
			#print("Visited node {:d} ({:d})".format(node.id, self._counter))
			#print(n, connections)
			if len(connections) > 0:
				#print("Adding {:d} connections".format(len(connections)))
				connected_ids.append([x.id for x in connections])
		sys.setrecursionlimit(lim)
		#print("recur out")
		return connected_ids

def predictedClusterSizeProbability(p, k):
	# From https://arxiv.org/pdf/1404.6806.pdf
	return (k**(k-2)/factorial(k)) * (2*p)**(k-1) * sp.exp(-2*k*p)


# Running into recursion errors
# Try networkx (https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.connected_components.html)

def findConnectedSubgraphs(node_indices, edge_indices):
	G = nx.Graph() # create graph
	nodeset = set(node_indices) # nodes in graph
	for (a,b) in edge_indices:
		G.add_edge(a,b)
	subgraphs = [[n for n in G.subgraph(c)] for c in nx.connected_components(G)] 
	# This will return only connected subgraphs; singletons will not be returned
	return subgraphs

