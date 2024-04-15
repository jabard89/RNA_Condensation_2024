#! python

import os, sys, math, unittest
# from [library] import [function]
import sedsim
import scipy as sp
import numpy as np

class test001(unittest.TestCase):
	def test_run(self):
		"""Initial setup"""
		molecules, mol_dict, molecule_ids = sedsim.generateMolecules(100)
		tube = sedsim.SedimentationVessel(10)
		tube.load(molecules)
		pos_before = [m.position_mm for m in tube.molecules]
		tube.spin(1e5, 10, 15*60)
		#for (i, mol) in enumerate(tube.molecules):
		#	print("{:d}\t{:1.2E}\t{:1.2E}\t{:1.2E}\t{:1.2E}\t{:1.2E}".format(
		#		mol.length, mol.mass, mol.radius, pos_before[i], mol.position, pos_before[i]-mol.position))

class test002(unittest.TestCase):
	def test_run(self):
		"""Terminal velocity"""
		vt = sedsim.getTerminalVelocity(mass=1e6/6.022e23, particle_radius=10e-9, rcf_g=100000, rotor_radius_cm=10, viscosity=1.002)
		self.assertTrue(vt > 8e-6 and vt < 9e-6)

class test003(unittest.TestCase):
	def test_run(self):
		"""Proportion pelleted"""
		molecules = []
		N = 20 # from 0 to N so actually N+1 types
		n_type = 1000
		min_length = 10
		max_length = 1e5
		labels = ["{:d}".format(n) for n in range(N+1)]
		mol_dict = {}
		for (n, nlab) in enumerate(labels):
			# log spacing
			L = np.exp(np.log(min_length) + (np.log(max_length)-np.log(min_length))*(n/N))
			# linear spacing
			#L = min_length + (max_length-min_length)*(n/N)
			radius = sedsim.getRNARadius(L)
			mol_dict[nlab] = []
			for nt in range(n_type):
				mol = sedsim.RNAMolecule(int(L), radius, nlab)
				molecules.append(mol)
				mol_dict[nlab].append(mol)
		tube = sedsim.SedimentationVessel(10)
		tube.load(molecules)
		#print("\nBefore:\n")
		for lab in labels:
			mols = mol_dict[lab]
			pelmols = [m for m in mols if m.pelleted]
			#n_label, n_pelleted, prop = tube.proportionPelleted(lab)
			#print("{:s}\t{:d}\t{:d}\t{:d}\t{:1.2f}".format(lab, mols[0].length, len(mols), len(pelmols), prop))
		tube.spin(1e5, 6, 15*60)
		#print("\nAfter:\n")
		for lab in labels:
			mols = mol_dict[lab]
			pelmols = [m for m in mols if m.pelleted]
			#n_label, n_pelleted, prop = tube.proportionPelleted(lab)
			#print("{:s}\t{:d}\t{:d}\t{:d}\t{:1.2f}".format(lab, mols[0].length, len(mols), len(pelmols), prop))

class test004(unittest.TestCase):
	def test_run(self):
		"""Clusters"""
		clust1 = sedsim.RNACluster()
		clust2 = sedsim.RNACluster()
		clust3 = sedsim.RNACluster()
		rna1 = sedsim.RNAMolecule(100, 1e-9)
		rna2 = sedsim.RNAMolecule(1000, 2e-9)
		rna3 = sedsim.RNAMolecule(10000, 3e-9)
		rna4 = sedsim.RNAMolecule(100000, 3e-9)
		clust1.add(rna1)
		clust1.add(rna2)
		clust2.add(rna3)
		clust2.add(clust1)
		clust3.add(clust2)
		self.assertTrue(clust3.count == 3)
		clust3.add(rna3)
		# rna3 already in clust3 via clust2, should not change number
		self.assertTrue(clust3.count == 3)
		# Adding new RNA should change numbers
		clust3.add(rna4)
		self.assertTrue(clust3.count == 4)
		self.assertTrue(clust2.count == 3)
		self.assertTrue(clust2.length == 11100)
		self.assertTrue(clust1.count == 2)
		self.assertTrue(clust1.length == 1100)

class test005(unittest.TestCase):
	def test_run(self):
		"""Clusters sedimenting"""
		clust1 = sedsim.RNACluster()
		clust2 = sedsim.RNACluster()
		rna0 = sedsim.RNAMolecule(10, 1e-9)
		rna1 = sedsim.RNAMolecule(100, 1e-9)
		rna2 = sedsim.RNAMolecule(1000, 2e-9)
		rna3 = sedsim.RNAMolecule(10000, 3e-9)
		clust1.add(rna1)
		clust1.add(rna2)
		clust2.add(rna3)
		clust2.add(clust1)
		molecules = [rna0, clust1, clust2]
		tube = sedsim.SedimentationVessel(10)
		tube.load(molecules)
		tube.spin(1e5, 6, 15*60)
		self.assertTrue(not rna0.pelleted)
		self.assertTrue(clust2.pelleted)

'''
class test006(unittest.TestCase):
	def test_run(self):
		"""Clusters merging"""
		rna0 = sedsim.RNAMolecule(10, 1e-9)
		rna1 = sedsim.RNAMolecule(100, 1e-9)
		rna2 = sedsim.RNAMolecule(1000, 2e-9)
		rna3 = sedsim.RNAMolecule(10000, 3e-9)
		col = sedsim.ClusterCollection()
		col.add(rna1)
		col.add(rna2)
		self.assertTrue(col.count==2)
		clust = col.merge(rna1, rna2)
		self.assertTrue(col.count==1)
		self.assertTrue(clust.length==1100)
		self.assertTrue(clust.radius>0.0)
'''

class test007(unittest.TestCase):
	def test_run(self):
		"""Crosslinking nucleotides"""
		rna = []
		col = sedsim.ClusterCollection()
		N = 10
		L = 1000
		for n in range(N):
			rnax = sedsim.RNAMolecule(L, 1, "{:d}".format(L))
			col.add(rnax)
		xlink = sedsim.DualCrosslinker(0.0, 1.0/L)
		self.assertTrue(col.count==N)
		new_col = xlink.crosslink(col)
		# should be clusters = fewer top-level items
		self.assertTrue(len([x for x in new_col]) < N)

class test008(unittest.TestCase):
	def test_run(self):
		"""Containment"""
		clust1 = sedsim.RNACluster()
		clust2 = sedsim.RNACluster()
		rna0 = sedsim.RNAMolecule(10, 1e-9)
		rna1 = sedsim.RNAMolecule(100, 1e-9)
		rna2 = sedsim.RNAMolecule(1000, 2e-9)
		rna3 = sedsim.RNAMolecule(10000, 3e-9)
		clust1.add(rna1)
		clust1.add(rna2)
		clust2.add(rna3)
		clust2.add(clust1)
		self.assertTrue(not rna3.isFree())
		self.assertTrue(clust2.isFree())
		self.assertTrue(not clust1.isFree())
		self.assertTrue(rna1.container == clust2)

class test009(unittest.TestCase):
	def test_containment(self):
		"""Containment and masses"""
		clust1 = sedsim.RNACluster()
		clust2 = sedsim.RNACluster()
		rna0 = sedsim.RNAMolecule(10, sedsim.getRNARadius(10))
		rna1 = sedsim.RNAMolecule(100, sedsim.getRNARadius(100))
		rna2 = sedsim.RNAMolecule(1000, sedsim.getRNARadius(1000))
		rna3 = sedsim.RNAMolecule(10000, sedsim.getRNARadius(10000))
		clust1.add(rna1)
		clust1.add(rna2)
		clust2.add(rna3)
		clust2.add(clust1)
		self.assertTrue(rna1.container == clust2)
		self.assertAlmostEqual(rna1.mass, sedsim.getRNAMass(10))
		self.assertAlmostEqual(clust1.mass, sedsim.getRNAMass(1100))
		self.assertAlmostEqual(clust2.mass, sedsim.getRNAMass(11100))

	def test_sedimentation(self):
		"""Containment and masses: sedimentation"""
		clust1 = sedsim.RNACluster()
		clust2 = sedsim.RNACluster()
		rna0 = sedsim.RNAMolecule(11100, sedsim.getRNARadius(11100))
		rna1 = sedsim.RNAMolecule(100, sedsim.getRNARadius(100))
		rna2 = sedsim.RNAMolecule(1000, sedsim.getRNARadius(1000))
		rna3 = sedsim.RNAMolecule(10000, sedsim.getRNARadius(10000))
		clust1.add(rna1)
		clust1.add(rna2)
		clust2.add(rna3)
		clust2.add(clust1)
		self.assertAlmostEqual(clust2.mass, sedsim.getRNAMass(11100))
		self.assertAlmostEqual(clust2.mass, rna0.mass)
		tube = sedsim.SedimentationVessel(10)
		tube.load([rna0, clust2])
		# Now sediment
		initial_pos = 10.0
		rna0.position_mm = initial_pos
		clust2.position_mm = initial_pos
		tube.spin(1e4, 6, 1*60)
		self.assertFalse(rna0.position_mm == initial_pos)
		self.assertFalse(clust2.position_mm == initial_pos)
		self.assertAlmostEqual(rna0.position_mm, clust2.position_mm)

class test010(unittest.TestCase):
	def test_crosslinking_per_molecule(self):
		"""Proportion free for per-molecule crosslinking"""
		mu = 0.3
		n = 10000
		col = sedsim.ClusterCollection()
		clink = sedsim.DualCrosslinker(mu,0)
		mols = [sedsim.RNAMolecule(100, sedsim.getRNARadius(100)) for i in range(n)]
		for mol in mols:
			col.add(mol)
		xcol = clink.crosslink(col)
		pfree = xcol.proportionFree()
		#print(pfree)
		self.assertAlmostEqual(pfree, 1.0-mu, 2)

	def test_crosslinking_per_length(self):
		"""Proportion free for per-nucleotide crosslinking"""
		nu = 0.01
		n = 1000
		L = 1000
		col = sedsim.ClusterCollection()
		clink = sedsim.DualCrosslinker(0,nu)
		mols = [sedsim.RNAMolecule(L, sedsim.getRNARadius(L)) for i in range(n)]
		for mol in mols:
			col.add(mol)
		xcol = clink.crosslink(col)
		pfree = xcol.proportionFree()
		#print(pfree)
		self.assertAlmostEqual(pfree, np.exp(-nu*L), 2)

class test011(unittest.TestCase):
	def test_crosslinking_per_molecule(self):
		"""Retrieve input parameters for per-molecule crosslinking"""
		px = None # probability of per-molecule crosslinking
		for xi in range(10):
			px = np.random.rand()
			n = 1000
			col = sedsim.ClusterCollection()
			clink = sedsim.DualVariableCrosslinker(px,0,1.0)
			mols = [sedsim.RNAMolecule(10, sedsim.getRNARadius(100)) for i in range(n)]
			for mol in mols:
				col.add(mol)
			xcol = clink.crosslink(col)
			pfree = xcol.proportionFree()
			sd_p = np.sqrt((1-px)*px/n) # binomial variance
			#print(pfree, 1-px, sd_p)
			# Is observed pfree within 2 SD of expectation = 1-px?
			self.assertTrue(pfree >= (1-px)-2*sd_p)
			self.assertTrue(pfree <= (1-px)+2*sd_p)
			#self.assertAlmostEqual(pfree, (1-px)*, 2)

class test012(unittest.TestCase):
	def test_psup_vs_pfree(self):
		px = np.random.rand()
		n = 1000
		col = sedsim.ClusterCollection()
		clink = sedsim.DualVariableCrosslinker(px,0,1.0)
		mols = [sedsim.RNAMolecule(10, sedsim.getRNARadius(100), 0.0) for i in range(n)]
		for mol in mols:
			col.add(mol)
		xcol = clink.crosslink(col)
		pfree = xcol.proportionFree()

class test013(unittest.TestCase):
	def test_molecule_conservation(self):
		np.random.seed(111)
		px = np.random.rand()
		n = 30
		col = sedsim.ClusterCollection()
		clink = sedsim.DualVariableCrosslinker(px,0,1.0)
		mols = [sedsim.RNAMolecule(10, sedsim.getRNARadius(100), 0.5) for i in range(n)]
		for mol in mols:
			col.add(mol)
		linked = clink.crosslink(col)
		#print(linked.count, col.count)
		self.assertTrue(set([m for m in linked.molecules]) == set([m for m in mols]))

class test014(unittest.TestCase):
	def test_particles(self):
		"""Full particles"""
		n = 4
		mols = [sedsim.RNAMolecule(10, sedsim.getRNARadius(100), 0.5) for i in range(n)]
		coll = sedsim.ClusterCollection()
		for m in mols:
			coll.add(m)
		coll_clust = sedsim.ClusterCollection()
		clust1 = sedsim.RNACluster(label='1')
		clust1.add(coll[0])
		clust1.add(coll[1])
		clust2 = sedsim.RNACluster(label='2')
		clust2.add(coll[2])
		clust2.add(coll[3])
		coll_clust.add(clust1)
		clust1.add(clust2)
		self.assertTrue(set([x for x in coll.molecules]) == set([x for x in coll_clust.molecules]))
		self.assertTrue(coll.count == 4)
		self.assertTrue(coll_clust.count == 4)
		# print("coll")
		# for (xi, m) in enumerate(coll):
		# 	print(xi)
		# print("coll_clust")
		# for (xi, m) in enumerate(coll_clust):
		# 	print(xi)

class test015(unittest.TestCase):
	def test_crosslinking_per_molecule(self):
		"""Proportion free for per-molecule crosslinking"""
		np.random.seed(111)
		mu = np.random.random()
		n = 4000
		col = sedsim.ClusterCollection()
		clink = sedsim.DualVariableCrosslinker(mu,0)
		mols = [sedsim.RNAMolecule(100, sedsim.getRNARadius(100), 0.0) for i in range(n)]
		mols = mols + [sedsim.RNAMolecule(100, sedsim.getRNARadius(100), 1.0) for i in range(n)]
		for mol in mols:
			col.add(mol)
		xcol = clink.crosslink(col)
		pfree = xcol.proportionFree()
		#print(pfree)
		self.assertAlmostEqual(pfree, 1.0-mu/2.0, 2)

class test015(unittest.TestCase):
	def test_loading_in_tube(self):
		"""Loading with tube"""
		n = 4
		mols = [sedsim.RNAMolecule(10, sedsim.getRNARadius(100), 0.5) for i in range(n)]
		coll = sedsim.ClusterCollection()
		for m in mols:
			coll.add(m)
		coll_clust = sedsim.ClusterCollection()
		clust1 = sedsim.RNACluster(label='1')
		clust1.add(coll[0])
		clust1.add(coll[1])
		clust2 = sedsim.RNACluster(label='2')
		clust2.add(coll[2])
		clust2.add(coll[3])
		coll_clust.add(clust1)
		clust1.add(clust2)
		# Make tube
		tube = sedsim.SedimentationVessel(20)
		tube.viscosity = 2
		tube.load(coll_clust)
		self.assertTrue(set([x for x in tube.molecules]) == set([x for x in coll.molecules]))
		self.assertTrue(coll.count == 4)
		self.assertTrue(tube.count == 4)

	def test_crosslinking_in_tube(self):
		"""Crosslinking with tube"""
		n = 1000
		L_short = 100
		L_long = 300000
		rcf_g = 2e4
		rotor_radius_cm = 20
		viscosity = 1.0
		spin_time_sec = 15*60
		tube_height_mm = 20
		coll = sedsim.ClusterCollection()
		clink = sedsim.DualVariableCrosslinker(0.5,0.0,1.0)
		mols = [sedsim.RNAMolecule(L_short, sedsim.getRNARadius(L_short), 0.0) for i in range(n)] + \
			[sedsim.RNAMolecule(L_long, sedsim.getRNARadius(L_long), 0.0) for i in range(n)]
		for mol in mols:
			coll.add(mol)
		xcol = clink.crosslink(coll)
		# Terminal velocity in m/s
		v1 = sedsim.getTerminalVelocity(mols[0].mass, mols[0].radius, rcf_g, rotor_radius_cm, viscosity)
		v2 = sedsim.getTerminalVelocity(mols[n].mass, mols[n].radius, rcf_g, rotor_radius_cm, viscosity)
		#print(tube_height_mm-v1*1000*spin_time_sec, \
		#	tube_height_mm-v2*1000*spin_time_sec)

		# Make tube
		tube = sedsim.SedimentationVessel(tube_height_mm)
		tube.viscosity = viscosity
		tube.load(xcol)
		tube.spin(rcf_g, rotor_radius_cm, spin_time_sec)
		self.assertTrue(set([x for x in tube.molecules]) == set([x for x in coll.molecules]))
		self.assertTrue(coll.count == 2*n)
		self.assertTrue(tube.count == 2*n)
		# Did anything spin down?
		pfree_short = 1-np.sum([m.pelleted for m in tube.molecules if m.length==L_short])/float(n)
		pfree_long = 1-np.sum([m.pelleted for m in tube.molecules if m.length==L_long])/float(n)
		#print(pfree_short, pfree_long)
		# Did the short things spin down less?
		self.assertTrue(pfree_short<1.0)
		self.assertTrue(pfree_long<pfree_short)

	def test_crosslinking_translation_per_molecule(self):
		"""Crosslinking with translation per molecule"""
		n = 100
		L_short = 10000
		rcf_g = 2e4
		rotor_radius_cm = 20
		viscosity = 1.0
		spin_time_sec = 15*60
		tube_height_mm = 20
		translation_prob = 1.0
		coll = sedsim.ClusterCollection()
		clink = sedsim.DualVariableCrosslinker(0.5,0.0,translation_prob)
		mols = [sedsim.RNAMolecule(L_short, sedsim.getRNARadius(L_short), 1.0) for i in range(n)] + \
			[sedsim.RNAMolecule(L_short, sedsim.getRNARadius(L_short), 0.0) for i in range(n)]
		for mol in mols:
			coll.add(mol)
		xcol = clink.crosslink(coll)

		# Make tube
		tube = sedsim.SedimentationVessel(tube_height_mm)
		tube.viscosity = viscosity
		tube.load(xcol)
		tube.spin(rcf_g, rotor_radius_cm, spin_time_sec)
		# Did anything spin down?
		pfree_trans = 1-np.sum([m.pelleted for m in mols[0:n]])/float(n)
		pfree_untrans = 1-np.sum([m.pelleted for m in mols[n:]])/float(n)
		#print(pfree_trans, pfree_untrans)
		# Did the short things spin down less?
		self.assertTrue(pfree_trans<1.0)
		self.assertTrue(pfree_untrans<pfree_trans)

	def test_crosslinking_translation(self):
		"""Crosslinking with translation"""
		n = 100
		L = 1000
		rcf_g = 2e4
		rotor_radius_cm = 20
		viscosity = 1.0
		spin_time_sec = 15*60
		tube_height_mm = 20
		full_translation = 1.0
		no_translation = 0.0
		per_molecule_translation = 0.9
		overall_pxlink_if_untranslated = 0.8

		# Full translation
		#print("#\n With translation ")
		coll1 = sedsim.ClusterCollection()
		clink1 = sedsim.DualVariableCrosslinker(overall_pxlink_if_untranslated,0.0, full_translation)
		mols1 = [sedsim.RNAMolecule(L, sedsim.getRNARadius(L), per_molecule_translation, "trans") for i in range(n)]
		for mol in mols1:
			coll1.add(mol)
		xcol1 = clink1.crosslink(coll1)

		# Zero translation
		#print("#\n No translation ")
		coll2 = sedsim.ClusterCollection()
		clink2 = sedsim.DualVariableCrosslinker(overall_pxlink_if_untranslated,0.0, no_translation)
		mols2 = [sedsim.RNAMolecule(L, sedsim.getRNARadius(L), per_molecule_translation, "untrans") for i in range(n)] 
		for mol in mols2:
			coll2.add(mol)
		xcol2 = clink2.crosslink(coll2)

		
		pfree_trans = np.mean([m.isFree() for m in mols1]) #/float(n)
		pfree_untrans = np.mean([m.isFree() for m in mols2]) #/float(n)
		print(pfree_trans, pfree_untrans)
		self.assertTrue(pfree_trans<1.0)
		# Did translated things cluster less?
		self.assertTrue(pfree_trans>pfree_untrans)

if __name__=="__main__":
	unittest.main(verbosity=2)
