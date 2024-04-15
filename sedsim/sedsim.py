import scipy as sp
import numpy as np
#from numpy.random import default_rng
import sys
import stats, util, na, dfstree
import collections, os

class SedimentableParticle(object):
	def __init__(self, label=None):
		self.label = label # type of protein, for example
		self.initial_position_mm = None
		self.position_mm = None
		self.pelleted = False
		self._psup = None # Proportion in the supernatant
		self._inside = None # pointer to cluster

	def addTo(self, containing_particle):
		res = False
		if self.container is not containing_particle.container:
			self._inside = containing_particle
			res = True
		return res

	def isFree(self):
		return self._inside is None

	def pellet(self):
		self.pelleted = True

	@property
	def psup(self):
		"""Proportion in the supernatant"""
		res = self.container._psup
		if res is None:
			# Update 
			res = 1.0
			if self.pelleted:
				res = 0.0
		return res

	@psup.setter
	def psup(self, psup_value):
		self._psup = psup_value

	@property
	def container(self):
		res = self
		if not self.isFree():
			res = self._inside
			while not res.isFree():
				res = res._inside
		return res

def getRNAMass(rna_length):
	# 330 g/mol nucleotide
	return rna_length * 330 / 6.022e23 # mass in grams

def getRNARadius(rna_length, scaling=0.34):
	# Source of this number
	return (2.96e-27*rna_length)**scaling

class RNAMolecule(SedimentableParticle):
	def __init__(self, length, radius, translation_prob=0.0, protection_prob=None, label=None):
		super().__init__(label)
		self.length = length
		self.radius = radius
		self.translation_prob = translation_prob
		self.protection_prob = protection_prob

	def isCompound(self):
		"""Is this a cluster of molecules"""
		return False

	@property
	def count(self):
		return 1
	
	@property
	def mass(self):
		return getRNAMass(self.length)

class RNACluster(SedimentableParticle):
	def __init__(self, label=None):
		super().__init__(label)
		self._particles = []
		self._total_length = 0
		self._numparticles = 0
		self._altered = False
		self._length = None
		self._radius = None
		self._mass = None
		self._count = None
		self._pelleted = False
		self._psup = None

	def add(self, particle):
		added = particle.addTo(self)
		if added:
			# if particle already belonged to this cluster, then don't do subsequent updates
			self._particles.append(particle)
			# DAD: use count? Or treat as lone particle?
			#self._numparticles += particle.count
			# Invalidate properties of this cluster so they will be recomputed
			self._length = None
			self._radius = None
			self._mass = None
			self._count = None

	# To avoid circular references where a cluster is added as a child
	# to a cluster that it contains, check lineages
	def lineage(self):
		"""Get parents back to top"""
		the_lineage = []
		if not self._parent is None:
			the_lineage = self._parent.lineage_recursive(the_lineage)
		the_lineage = [self] + the_lineage
		return the_lineage

	def lineage_recursive(self, the_lineage):
		"""Turn hierarchy into a list"""
		the_lineage.append(self)
		if not self._parent is None:
			the_lineage += self._parent.lineage_recursive(the_lineage)
		return the_lineage

	def pellet(self):
		super().pellet()
		for m in self:
			m.pellet()

	@property
	def count(self):
		"""The number of molecules in the cluster"""
		if self._count is None:
			self._count = len([x for x in self.molecules])
		return self._count
		#return len([x for x in self.particles])

	@property
	def length(self):
		"""The total length of the cluster"""
		if self._length is None:
			self._length = sum([m.length for m in self.molecules])
		return self._length

	@property
	def mass(self):
		if self._mass is None:
			self._mass = sum([m.mass for m in self.molecules])
		return self._mass

	@property
	def radius(self):
		if self._radius is None:
			self._radius = getRNARadius(self.length)
		return self._radius

	def isCompound(self):
		"""Is this a cluster of molecules"""
		return True

	@property
	def molecules(self):
		"""Traverse the tree structure"""
		for m in self._particles:
			if m.isCompound():
				for mchild in m.molecules:
					yield mchild
			else:
				yield m

	def __iter__(self):
		for m in self.molecules:
			yield m

class ClusterCollection(object):
	"""A collection of clusters. 
	"""
	def __init__(self):
		self._particles = []
		self._count = None

	def add(self, particle):
		self._particles.append(particle)
		# Invalidate count so that it will be recomputed
		self._count = None

	def remove(self, particle):
		self._particles.remove(particle)
		# Invalidate count so that it will be recomputed
		self._count = None

	# Pelleting
	def pellet(self):
		for p in self:
			p.pellet()

	# Count collection
	@property
	def count(self):
		"""The number of clusters in the collection."""
		return len(self._particles)

	@property
	def count(self):
		"""The number of molecules in the cluster"""
		if self._count is None:
			self._count = len([x for x in self.molecules])
		return self._count

	# Count collection
	def __len__(self):
		"""The number of particles in the collection."""
		return self.count

	# # Iterate over the top-level collection
	# def items(self):
	# 	for p in self:
	# 		yield p

	# Iterate depth-first over full collection
	@property
	def molecules(self):
		for p in self._particles:
			if p.isCompound():
				for p2 in p.molecules:
					yield p2
			else:
				yield p

	def __iter__(self):
		for p in self._particles:
			yield p

	def proportionFree(self):
		tot = 0.0
		free = 0.0
		for m in self.molecules:
			tot += 1
			if m.isFree():
				free += 1
		return free/tot


	def __getitem__(self, index):
		assert(index < len(self._particles))
		return self._particles[index]

class SedimentationVessel(ClusterCollection):
	"""
	A centrifuge tube by another name.
	"""

	def __init__(self, height_mm):
		super().__init__()
		self.height_mm = height_mm

		# viscosity
		self._viscosity = 1.002 # g /(m . sec), default 1.002 for water
	
	@property
	def viscosity(self):
		return self._viscosity

	@viscosity.setter
	def viscosity(self, eta):
		self._viscosity = eta

	def load(self, collection):
		self._particles = [x for x in collection]
		for part in self._particles:
			# Choose uniform distribution over tube
			part.position_mm = np.random.random()*self.height_mm
			part.initial_position_mm = part.position_mm

	def spin(self, rcf_g, rotor_radius_cm, time):
		"""Do physics"""
		# For mo
		# vt = (m r omega^2)/(6 pi eta r_0)
		# accel = (r omega^2)
		# vt = (m/r_0) * accel/(6 pi eta)
		# 
		# g force (RCF) = 1.118e-5 * RPM^2 * r_cm
		# RPM = sqrt(RCF * 89445.44 / r_cm)
		# omega = radians/sec = rev / minute * (1 min/60 sec) * (2 pi radians/rev)
		#
		# Typical values
		# r = 0.1m
		# rpm = 30000; omega = 3141
		# m = 1 MDa = 1e6 g/mol * (1 mol/ 6.022e23 molecules) = 1.66e-18 g
		# r0 = 10 nm
		# eta = 1.002 g/(m.sec) for water
		# vt = (1.66e-18 / 10e-9) * 0.1 * 3141^2 / (6 * pi * 1.002)
		#    = 8.67 e -6 m/sec = 8.67 e -3 mm/sec =~ 0.01 mm/sec so will move one mm in 100 sec, not unreasonable
		
		# convert RCF (force) to RPM (angular velocity)
		rpm = RPMfromRCFg(rcf_g, rotor_radius_cm)
		omega = rpm * (2 * sp.pi)/60 # rotations per minute --> radians per second
		accel = 0.01 * rotor_radius_cm * omega*omega # m /sec^2
		prefactor = accel/(6 * sp.pi * self._viscosity) # (m/s^2) / (g/(m . s)) = m^2/(g . s)
		# Iterate over particles
		for m in self:
			terminal_vel = prefactor * m.mass/m.radius # meters/sec
			#print(m.mass, terminal_vel, vt)
			diffusion = 0.0 # DAD: add diffusion
			if not m.pelleted:
				# Velocity is in meters/sec, move to mm/sec
				# Multiply by time to get distance
				new_pos = m.position_mm - terminal_vel*1000*time + diffusion
				if new_pos < 0.0: # Particle hit the bottom of the tube
					new_pos = 0.0
					m.pellet()
					#m.pelleted = True
				# Update position
				m.position_mm = new_pos

class DeterministicSedimentationVessel(SedimentationVessel):
	"""
	A centrifuge tube by another name.
	In this tube, particles obey average behavior and there is no
	random component to their sedimentation
	"""

	def __init__(self, height_mm):
		super().__init__(height_mm)

	def load(self, collection):
		self._particles = sorted([x for x in collection], key=lambda x: x.length)

		#for part in self._particles:
		#	# Choose uniform distribution over tube
		#	part.position_mm = np.random.random()*self.height_mm
		#	part.initial_position_mm = part.position_mm

	def spin(self, rcf_g, rotor_radius_cm, time):
		"""Do physics"""
		# For mo
		# vt = (m r omega^2)/(6 pi eta r_0)
		# accel = (r omega^2)
		# vt = (m/r_0) * accel/(6 pi eta)
		# 
		# g force (RCF) = 1.118e-5 * RPM^2 * r_cm
		# RPM = sqrt(RCF * 89445.44 / r_cm)
		# omega = radians/sec = rev / minute * (1 min/60 sec) * (2 pi radians/rev)
		#
		# Typical values
		# r = 0.1m
		# rpm = 30000; omega = 3141
		# m = 1 MDa = 1e6 g/mol * (1 mol/ 6.022e23 molecules) = 1.66e-18 g
		# r0 = 10 nm
		# eta = 1.002 g/(m.sec) for water
		# vt = (1.66e-18 / 10e-9) * 0.1 * 3141^2 / (6 * pi * 1.002)
		#    = 8.67 e -6 m/sec = 8.67 e -3 mm/sec =~ 0.01 mm/sec so will move one mm in 100 sec, not unreasonable
		
		# convert RCF (force) to RPM (angular velocity)
		rpm = RPMfromRCFg(rcf_g, rotor_radius_cm)
		omega = rpm * (2 * sp.pi)/60 # rotations per minute --> radians per second
		accel = 0.01 * rotor_radius_cm * omega*omega # m /sec^2
		prefactor = accel/(6 * sp.pi * self._viscosity) # (m/s^2) / (g/(m . s)) = m^2/(g . s)
		# Iterate over particles
		for m in self:
			terminal_vel = prefactor * m.mass/m.radius # meters/sec
			#print(m.mass, terminal_vel, vt)
			dist_moved = terminal_vel*1000*time
			# For a uniformly distributed particle in a tube of height h, and distance moved d at terminal
			# velocity, a fraction d/h will pellet.
			frac_pelleted = min(dist_moved/self.height_mm, 1.0)
			m.psup = 1.0-frac_pelleted


###
# 
# Crosslinkers
#
####

class DualCrosslinker(object):
	"""Crosslink molecules and nucleotides"""
	def __init__(self, mol_crosslink_prob, nt_crosslink_prob):
		# Crosslink ratio = ratio of markers (crosslink ends) to molecules, which can be >1
		self.mol_crosslink_prob = mol_crosslink_prob
		self.nt_crosslink_prob = nt_crosslink_prob
		#print("P = ({:f},{:f})".format(self.mol_crosslink_prob, self.nt_crosslink_prob))

	# @util.printTiming
	def crosslink(self, cluster_collection):
		"""Crosslink molecules with a given probability per nucleotide and average ."""
		# Count nucleotides, N
		coll = cluster_collection
		num_molecules = coll.count
		if num_molecules < 2:
			return coll

		# Nothing to do if probability of crosslinking is zero.
		if self.mol_crosslink_prob==0.0 and self.nt_crosslink_prob==0.0:
			return coll

		#print("Now crosslinking")

		crosslink_pairs = []

		# Molecule crosslinking
		# Generate crosslinks such that the probability that a molecule is crosslinked is mol_crosslink_prob (p)
		# Balls and buckets: throw k balls into N buckets and the fraction filled is p = 1-(1-1/N)^k
		# so must throw k = log(1-p)/log(1-1/N) balls to fill a fraction p on average
		# In this case each "ball" is a crosslink which hits two "buckets" -- molecules -- so that
		# to get p molecules crosslinked we through k/2 crosslinks
		# 1) Go through list of molecules and select each with probability p into a crosslink list
		# 2) For 
		# 1) Given N molecules, sample M ~ N*p
		# 2) For each molecule m_i in the list, link it to some random molecule m_j in the list.
		# 3) Done. Now (on average) N*p molecules have a crosslink, and clusters can form.

		# Nucleotide crosslinking
		# Generate crosslinks such that the probability that a nucleotide is crosslinked is nuc_crosslink_prob (p)
		# 1) Given N molecules with lengths l_i


		# Generate a total number of crosslinks such that
		# the ratio of crosslinks to molecules is = mol_crosslink_prob
		num_crosslinks = int(0.5*np.ceil(np.log(1-self.mol_crosslink_prob)/np.log(1-1.0/num_molecules)))
		# Choose molecules for each crosslink

		#num_molecules_to_crosslink = int(np.ceil(self.mol_crosslink_prob * num_molecules))
		#print("Generating {:d} crosslinks".format(num_crosslinks))
		# Choose sequences to crosslink
		# Self-crosslinks not allowed
		molecule_indices = [i for i in range(num_molecules)]
		rng = np.random.default_rng()
		mol_pair_list = [rng.choice(molecule_indices, 2, replace=False) for xi in range(num_crosslinks)]
		#print("Generating {:d} mol crosslinks (prob = {:E})".format(num_crosslinks, self.mol_crosslink_prob))
		crosslink_pairs += mol_pair_list
		#print(len(crosslink_pairs))

		# Nucleotide crosslinking
		crosslink_coll = coll
		if self.nt_crosslink_prob>0.0:
			lengths = [m.length for m in crosslink_coll]
			# Total number of nucleotides
			num_nucleotides = int(sum(lengths))
			# Generate a total number of crosslinks C such that
			# the ratio of nucleotide crosslinks to nucleotides is = nt_crosslink_prob
			num_nt_crosslinks = int(np.ceil(num_nucleotides * self.nt_crosslink_prob/2.0))
			#print("Generating {:d} nt crosslinks (prob = {:E})".format(num_nt_crosslinks, self.nt_crosslink_prob))
			#print("Generating {:d} crosslinks".format(num_nt_crosslinks))
			# Choose nucleotides to crosslink 
			# Meaning choose sequences, weighted by length
			norm_weights = [L/num_nucleotides for L in lengths]
			# Choose nucleotides to crosslink = choose sequences by length
			# We allow self-crosslinking here
			ntx_seq_list = list(rng.choice(molecule_indices, num_nt_crosslinks, replace=True, p=norm_weights))
			# Turn nts into edges
			half_n = int(np.floor(len(ntx_seq_list)/2))-1
			ntx_pair_list = [(ntx_seq_list[2*i], ntx_seq_list[2*i+1]) for i in range(half_n)]
			# Add nucleotide crosslinks to list
			crosslink_pairs += ntx_pair_list

		# Find connected clusters
		connected_sets = dfstree.findConnectedSubgraphs(molecule_indices, crosslink_pairs)
		unconnected_seqs = set(molecule_indices).difference(set([x for xs in connected_sets for x in xs]))

		# rebuild cluster collection
		new_cluster_collection = ClusterCollection()
		for (setind, conset) in enumerate(connected_sets):
			clust = RNACluster(label='{:d}'.format(setind))
			for id in conset:
				mol = coll[id]
				clust.add(mol)
			new_cluster_collection.add(clust)
		for seqind in unconnected_seqs:
			# singletons, add directly
			mol = coll[seqind]
			new_cluster_collection.add(mol)
		return new_cluster_collection

class DualVariableCrosslinker(object):
	"""Crosslink molecules and nucleotides with variation based on translation status """
	def __init__(self, mol_crosslink_prob, nt_crosslink_prob, translation_prob):
		# Crosslink ratio = ratio of markers (crosslink ends) to molecules, which can be >1
		self.mol_crosslink_prob = mol_crosslink_prob
		# Calculate a Poisson rate mu: P(xlink) = 1-P(0 xlink) = 1-exp(-mu)
		# so mu = -log(1-P(xlink))
		self.mol_crosslink_rate = -np.log(1-self.mol_crosslink_prob)
		self.nt_crosslink_prob = nt_crosslink_prob
		# Similarly compute rate for crosslinking nucleotides
		self.nt_crosslink_rate = -np.log(1-self.nt_crosslink_prob)
		self.translation_prob = translation_prob
		#print("P = ({:f},{:f})".format(self.mol_crosslink_prob, self.nt_crosslink_prob))

	def crosslink(self, cluster_collection):
		"""Crosslink molecules with a given probability per nucleotide and average ."""
		# Count nucleotides, N
		coll = cluster_collection
		total_molecules = coll.count
		if total_molecules < 2:
			# Can't crosslink 1 or fewer molecules
			return coll

		# Nothing to do if probability of crosslinking is zero.
		if self.mol_crosslink_prob==0.0 and self.nt_crosslink_prob==0.0:
			return coll

		# Identify molecules to be crosslinked. Here we will select molecules based on the probability
		# that they have zero ribosomes. Each molecule has a translation probability, and there is 
		# an overall translation probability for the cell which multiplies the molecule-specific probability.
		# Note that for this model, full translation (both values = 1.0) prevents any crosslinking.
		molecule_indices = []
		crosslink_coll = ClusterCollection()
		for xi,mol in enumerate(coll):
			if np.random.random() > self.translation_prob*mol.translation_prob*(1.0-mol.protection_prob): # will hit with probability 1-self.trans_prob*mol.trans_prob*(1-mol.protection_prob)
				molecule_indices.append(xi)
				crosslink_coll.add(mol)
		num_molecules = len(molecule_indices)

		if num_molecules == 0:
			# No molecules competent for crosslinking. Just return original list
			return coll

		crosslink_pairs = []

		# Molecule crosslinking
		# Generate crosslinks such that the probability that a molecule is crosslinked is mol_crosslink_prob (p)
		# Balls and buckets: throw k balls into N buckets and the fraction filled is p = 1-(1-1/N)^k
		# so must throw k = log(1-p)/log(1-1/N) balls to fill a fraction p on average
		# In this case each "ball" is a crosslink which hits two "buckets" -- molecules -- so that
		# to get p molecules crosslinked we through k/2 crosslinks
		# DAD: special-case self.mol_crosslink_prob = 1? To do.
		num_crosslinks = int(0.5*np.ceil(np.log(1-self.mol_crosslink_prob)/np.log(1-1.0/num_molecules)))
		# Choose molecules for each crosslink

		# Choose sequences to crosslink
		# Self-crosslinks not allowed
		#molecule_indices = [i for i in range(num_molecules)]
		rng = np.random.default_rng()
		mol_pair_list = [tuple(rng.choice(molecule_indices, 2, replace=False)) for xi in range(num_crosslinks)]
		#print("Generating {:d} mol crosslinks (prob = {:E})".format(num_crosslinks, self.mol_crosslink_prob))
		crosslink_pairs += mol_pair_list
		#print(len(crosslink_pairs))
		#print(mol_pair_list)

		# Nucleotide crosslinking
		if self.nt_crosslink_prob>0.0:
			lengths = [m.length for m in crosslink_coll]
			# Total number of nucleotides
			num_nucleotides = int(sum(lengths))
			# Generate a total number of crosslinks C such that
			# the ratio of nucleotide crosslinks to nucleotides is = nt_crosslink_prob
			num_nt_crosslinks = int(np.ceil(num_nucleotides * self.nt_crosslink_prob/2.0))
			#print("Generating {:d} nt crosslinks (prob = {:E})".format(num_nt_crosslinks, self.nt_crosslink_prob))
			#print("Generating {:d} crosslinks".format(num_nt_crosslinks))
			# Choose nucleotides to crosslink 
			# Meaning choose sequences, weighted by length
			norm_weights = [L/num_nucleotides for L in lengths]
			# Choose nucleotides to crosslink = choose sequences by length
			# We allow self-crosslinking here
			ntx_seq_list = list(rng.choice(molecule_indices, num_nt_crosslinks, replace=True, p=norm_weights))
			# Turn nts into edges
			half_n = int(np.floor(len(ntx_seq_list)/2))-1
			# print(len(ntx_seq_list), half_n)
			# ntx_pair_list = []
			# for i in range(half_n):
			# 	nsi = ntx_seq_list[2*i]
			# 	nsj = ntx_seq_list[2*i+1]
			# 	#mi = molecule_indices[nsi]
			# 	#mj = molecule_indices[nsj]
			# 	ntx_pair_list.append((mi, mj))
			ntx_pair_list = [(ntx_seq_list[2*i], ntx_seq_list[2*i+1]) for i in range(half_n)]
			# Add nucleotide crosslinks to list
			crosslink_pairs += ntx_pair_list

		# Find connected clusters
		connected_sets = dfstree.findConnectedSubgraphs(molecule_indices, crosslink_pairs)
		unconnected_seqs = set([x for x in range(total_molecules)]).difference(set([x for xs in connected_sets for x in xs]))

		# rebuild cluster collection
		new_cluster_collection = ClusterCollection()
		for (setind, conset) in enumerate(connected_sets):
			clust = RNACluster(label='{:d}'.format(setind))
			for seqind in conset:
				clust.add(coll[seqind])
			new_cluster_collection.add(clust)
		for seqind in unconnected_seqs:
			# singletons, add directly
			new_cluster_collection.add(coll[seqind])
		return new_cluster_collection


# Summary setup
class MoleculeSummary(object):
	"""Store average information about molecules"""
	def __init__(self):
		self.mass = None
		self.length = None
		# Proportion in supernatant, not pellet
		self.psup = None
		# Proportion free (not in pellet and not in cluster)
		self.pfree = None
		# Proportion clustered (whether or not clusters are pelleted)
		self.pclust = None
		self.id = None
		self.count = None
		self.cluster_count = None
		self.cluster_mass = None
		self.cluster_radius = None


####
#
#  Utility functions
#
####

class Template(object):
	"""A place to store transcript identifier, abundance, length, and other features."""
	def __init__(object):
		id = None
		abund = None
		length = None
		translation = None
		viscosity = None
		rcf_g = None # Rotor centrifugal force, in g


def RPMfromRCFg(rcf_g, rotor_radius_cm):
	# g force (RCF) = 1.118e-5 * RPM^2 * r_cm
	# RPM = sqrt(RCF * 89445.44 / r_cm)
	rpm = np.lib.scimath.sqrt(rcf_g * 89445.44 / rotor_radius_cm)
	return rpm

def getTerminalVelocity(mass, particle_radius, rcf_g, rotor_radius_cm, viscosity):
	# vt = (m r omega^2)/(6 pi eta r_0)
	rpm = RPMfromRCFg(rcf_g, rotor_radius_cm)
	omega = rpm * (2 * sp.pi)/60
	vt = (mass/particle_radius) * rotor_radius_cm*0.01 * omega*omega/(6 * sp.pi * viscosity)
	return vt

class ControlSettings(object):
	def __init__(self):	
		self.number_of_molecules = None
		self.crosslink_method = None
		self.experiment_id = None
		self.template_id_field = None
		self.template_filename = None
		self.template_length_field = None
		self.template_abundance_field = None
		self.template_translation_field = None
		self.template_protection_field = None
		self.crosslink_mol_prob = None
		self.crosslink_nt_prob = None
		self.translation_prob = None
		self.viscosity_mpa_sec = None
		self.spin_time_sec = None
		self.sample_height_mm = None
		self.rotor_radius_cm = None
		self.rcf_g = None
		self.replicates = None
		self.random_seed = None

def loadControlFile(inf):
	"""Load control file"""
	dlr = util.DelimitedLineReader(inf, header=True)
	settings_dict = {}
	for flds in dlr.dictentries:
		settings = ControlSettings()
		settings.number_of_molecules = flds['num.molecules']
		settings.crosslink_method = flds['crosslink.method']
		settings.experiment_id = flds['experiment.id']
		settings.template_id_field = flds['template.id.field']
		settings.template_filename = flds['template.filename']
		settings.template_length_field = flds['template.length.field']
		settings.template_abundance_field = flds['template.abundance.field']
		settings.template_translation_field = flds['template.translation.field']
		settings.template_protection_filename = flds['template.protection.filename']
		settings.template_protection_field = flds['template.protection.field']
		settings.crosslink_mol_prob = flds['crosslink.mol.prob']
		settings.crosslink_nt_prob = flds['crosslink.nt.prob']
		settings.translation_prob = flds['translation.prob']
		settings.viscosity_mpa_sec = flds['viscosity.mpa.sec']
		settings.spin_time_sec = flds['spin.time.sec']
		settings.sample_height_mm = flds['sample.height.mm']
		settings.rotor_radius_cm = flds['rotor.radius.cm']
		settings.rcf_g = flds['rcf.g']
		settings.replicates = flds['replicates']
		settings.random_seed = flds['random.seed']
		settings_dict[settings.experiment_id] = settings
	return settings_dict


def generateMoleculesFromTemplates(template_settings):
	mol_dict = {}
	molecules = []
	templates = {}
	# Background data keyed by template_id_field (e.g. ORF)
	lengths = {}
	abundances = {}
	translation_probs = {}
	protection_probs = {}


	# DAD: for now, assume all the settings are correct.
	lengths = readTemplateFile(template_settings.template_filename, template_settings.template_id_field, template_settings.template_length_field)
	abundances = readTemplateFile(template_settings.template_filename, template_settings.template_id_field, template_settings.template_abundance_field)
	#translation_probs = readTemplateFile(template_settings.template_filename, template_settings.template_id_field, template_settings.template_translation_field)
	#protection_probs = readTemplateFile(template_settings.template_protection_filename, template_settings.template_id_field, template_settings.template_protection_field)

	template_fname = os.path.expanduser(template_settings.template_filename)
	# Read input
	if not os.path.isfile(template_fname):
		raise IOError("# Error: file {} does not exist".format(template_fname))
	with open(template_fname,'r') as inf:
		# Read a tab-delimited file?
		dlr = util.DelimitedLineReader(inf, header=True)
		#idfld = template_settings.template_id_field
		#abundfld = template_abundance_field
		#lengthfld = template_length_field
		#tlnfld = template_translation_field

		# Iterate over genes
		for flds in dlr.dictentries:
			template = Template()
			template.id = flds[template_settings.template_id_field]
			template.length = lengths[template.id]
			template.abund = abundances[template.id]
			template.translation = 1.0
			template.protection = 0
			#try:
			#	template.translation = translation_probs[template.id]
			#except KeyError:
			#	template.translation = None
			#if na.isNA(template.translation):
			#	template.translation = 0.85 # DAD: should be calculated from data; this is median of 30C occupancy
			#if template.translation < 0.0 or template.translation > 1.0:
			#	raise ValueError("# Error: translation probability outside (0,1) interval (= {})".format(template.translation))
			#try:
			#	template.protection = protection_probs[template.id]
			#except KeyError:
			#	template.protection = None
			# template.abund = flds[abundfld]
			# template.translation = 1.0 # Default to full translation
			# if not tlnfld is None:
			# 	template.translation = flds[tlnfld]
			# 	if na.isNA(template.translation):
			# 		template.translation = 0.85 # DAD: should be calculated from data; this is median of 30C occupancy
			# 	if template.translation < 0.0 or template.translation > 1.0:
			# 		raise ValueError("# Error: translation probability outside (0,1) interval (= {})".format(template.translation))
			# if not tlnfld is None:
			# 	template.translation = flds[tlnfld]
			# 	if na.isNA(template.translation):
			# 		template.translation = 0.85 # DAD: should be calculated from data; this is median of 30C occupancy
			# 	if template.translation < 0.0 or template.translation > 1.0:
			# 		raise ValueError("# Error: translation probability outside (0,1) interval (= {})".format(template.translation))
			# if na.isNA(template.abund): # What to do if abundance is NA? Set to 1/number of molecules.
			# 	template.abund = 1.0/number_of_molecules
			# template.length = flds[lengthfld]
			if not na.isNA(template.abund) and not na.isNA(template.length):
				templates[template.id] = template
	# Strategy: don't attempt to precisely match abundance (due to fractional transcripts per cell etc.).
	# Instead, sample N transcripts according to the abundance distribution.
	# DAD: may want to require that certain transcripts are present in at least one copy
	total_abundance = sum([t.abund for (id,t) in templates.items()])
	molecule_ids = []
	for (id,t) in templates.items():
		molecule_ids.append(id)
		p_abund = t.abund/total_abundance
		#print("{},{:e},{},{}".format(id, p_abund, t.abund, total_abundance))
		num_to_generate = int(p_abund * template_settings.number_of_molecules)
		if num_to_generate > 0:
			mol_dict[id] = []
			for n in range(num_to_generate):
				# DAD: suppress translation and protection in v1
				mol = RNAMolecule(int(t.length), getRNARadius(t.length), t.translation, t.protection, t.id)
				#mol = RNAMolecule(int(t.length), getRNARadius(t.length), None, None, t.id)
				molecules.append(mol)
				mol_dict[id].append(mol)
	return molecules, mol_dict, molecule_ids

def readTemplateFile(filename, template_id_field, template_value_field):
	template_fname = os.path.expanduser(filename)
	# Read input
	value_dict = {}
	if not os.path.isfile(template_fname):
	 	raise IOError("# Error: file {} does not exist".format(template_fname))
	with open(template_fname,'r') as inf:
		dlr = util.DelimitedLineReader(inf, header=True)
		for flds in dlr.dictentries:
			id = flds[template_id_field]
			value = flds[template_value_field]
			# DAD: should handle cases when value is NA
			value_dict[id] = value
	return value_dict


def generateMolecules(n_types, n_type=10, min_length=10, max_length=1e5):
	# Generate de novo
	mol_dict = {}
	molecules = []
	N = n_types
	labels = ["{:d}".format(n) for n in range(N+1)]
	molecule_ids = labels
	for (n, nlab) in enumerate(labels):
		# log spacing
		L = np.exp(np.log(min_length) + (np.log(max_length)-np.log(min_length))*(n/N))
		# linear spacing
		#L = min_length + (max_length-min_length)*(n/N)
		radius = getRNARadius(L)
		translation = 1.0 # probability of translation
		mol_dict[nlab] = []
		for nt in range(n_type):
			mol = RNAMolecule(int(L), radius, translation, nlab)
			#print(mol.length, mol.radius, mol.mass)
			molecules.append(mol)
			mol_dict[nlab].append(mol)
	return molecules, mol_dict, molecule_ids
