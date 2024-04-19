#! python

import sys, os, math, random, argparse
#print(sys.version)
#print(sys.path)
import util, biofile, na
from collections import defaultdict
import scipy as sp
import numpy as np
import sedsim

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Simulation of sedimentation")
	# Required arguments
	#parser.add_argument(dest="in_fname", default=None, help="input filename")
	# Optional arguments
	parser.add_argument("--experiment-file", dest="experiment_file", type=str, default=None, help="experiment file for parameters")
	# The following arguments can override values in the experiment file provided above.
	#parser.add_argument("--template-file", dest="template_fname", default=None, help="molecular template filename, tab-delimited")
	#parser.add_argument("--template-abundance-field", dest="template_abundance_fld", default=None, help="column in template file for molecular abundance")
	#parser.add_argument("--template-id-field", dest="template_id_fld", default=None, help="column in template file for molecule species identifier")
	#parser.add_argument("--template-length-field", dest="template_length_fld", default=None, help="column in template file for length of molecule")
	#parser.add_argument("--template-translation-field", dest="template_translation_fld", default=None, help="column in template file for translation of molecule")
	#parser.add_argument("--template-protection-field", dest="template_protection_fld", default=None, help="column in template file for protection of molecule")

	# Spin conditions
	parser.add_argument("--rcf-g", dest="rcf_g", type=float, default=20000, help="rotor centrifugal force in g")
	parser.add_argument("--rotor-radius", dest="rotor_radius_cm", type=float, default=10, help="rotor radius in centimeters")
	parser.add_argument("--sample-height", dest="sample_height_mm", type=float, default=10, help="sample height in millimeters")
	parser.add_argument("--spin-time", dest="spin_time_sec", type=float, default=15*60, help="spin time in seconds")
	parser.add_argument("--viscosity", dest="viscosity", type=float, default=1.002, help="viscosity in g/m.sec")

	parser.add_argument("--number-of-molecules", dest="number_of_molecules", type=int, default=1000, help="number of molecules")

	parser.add_argument("--reps", dest="replicates", type=int, default=1, help="number of replicates to generate")
	parser.add_argument("--seed", dest="random_seed", type=int, default=111, help="random number generator seed")

	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("--summary-only", dest="summary_only", action="store_true", default=False, help="don't produce full output file")
	parser.add_argument("--summary-out", dest="summary_out_fname", default=None, help="summary-output filename")

	#parser.add_argument("--viscosity", dest="viscosity", default=1.002, help="viscosity")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	header_outs = util.OutStreams()
	if not options.summary_only:
		data_outs = util.OutStreams()
		header_outs.addStream(data_outs)
	summary_outs = util.OutStreams()
	
	# Start up output
	if not options.summary_only:
		if not options.out_fname is None:
			outf = open(options.out_fname,'w')
			data_outs.addStream(outf)
		else:
			# By default, write to stdout
			data_outs.addStream(sys.stdout)

	summary_output = not options.summary_out_fname is None
	if summary_output:
		summary_outf = open(options.summary_out_fname,'w')
		summary_outs.addStream(summary_outf)
		header_outs.addStream(summary_outs)

	# Write out parameters
	header_outs.write("# Run started {}\n".format(util.timestamp()))
	header_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	header_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		header_outs.write("#\t{k}: {v}\n".format(k=k, v=v))


	# Load control file
	# DAD: load with sedsim.loadControlFile()
	# Set seed etc. from these parameters

	experiment_settings_dict = {}
	experiment_order = []
	if not options.experiment_file is None:
		exp_filename = os.path.expanduser(options.experiment_file)
		if not os.path.isfile(exp_filename):
		 	raise IOError("# Error: file {} does not exist".format(exp_filename))
		with open(exp_filename,'r') as inf:
		 	experiment_settings_dict = sedsim.loadControlFile(inf)
		 	# DAD: implement ordering by experiment number
		 	experiment_order = [k for k in experiment_settings_dict]
	#print(experiment_order)

	# Essential idea: simulate sedimentation and crosslinking.
	# Load tube with mRNAs, bound by protein, crosslinked
	# Spin tube to terminal velocity
	# Count mRNAs in pellet and in supernatant

	#  python sedsim-driver.py --crosslink-prob 1e-5 --template-file ./scer-orf-tpm-occ-length.txt --template-id-field ORF --template-abundance-field mrna --template-length-field LengthTxEst --template-translation-field Occ.mean --out test.txt

	# Write output
	dout = util.DelimitedOutput()
	dout.addHeader('experiment.id','experiment identifier','s')
	dout.addHeader('spin.rcf.g','rotor centrifugal force in g','e')
	dout.addHeader('spin.time.sec','ratio of crosslinks to molecules','f')
	dout.addHeader('replicate','replicate number','d')
	dout.addHeader('crosslink.mol.prob','ratio of crosslinks to molecules','f')
	dout.addHeader('crosslink.nt.prob','ratio of crosslinks to nucleotides','e')
	dout.addHeader('translation.prob','fraction of maximal translation rate','e')
	dout.addHeader('mol.translation.prob','probability this molecule is translated at maximal translation rate','e')
	dout.addHeader('mol.protection.prob','probability this molecule is protected from crosslinking','e')
	dout.addHeader('id','molecule ID','s')
	dout.addHeader('length','length in nucleotides','d')
	dout.addHeader('mass','mass in g/mol','f')
	dout.addHeader('radius','radius in nanometers','f')
	dout.addHeader('psup','average proportion in the supernatant','f')
	dout.addHeader('cluster.id','unique identifier of cluster','s')
	dout.addHeader('cluster.size','number of molecules in cluster','d')
	dout.addHeader('cluster.length','cluster length in nucleotides','d')
	dout.addHeader('cluster.mass','cluster mass in g/mol','e')
	dout.addHeader('cluster.radius','cluster radius in nanometers','f')
	dout.addHeader('cluster.initial.position','cluster initial position in mm above tube bottom','f')
	dout.addHeader('cluster.position','cluster position in mm above tube bottom','f')
	dout.addHeader('num.pelleted','number of molecules hit the bottom of the tube','d')
	dout.addHeader('num.total','number of molecules','d')
	#dout.addHeader('prop.pelleted.noiseless','proportion of molecules that hit the bottom of the tube','f')
	#dout.addHeader('prop.pelleted','proportion of molecules at the bottom of the tube','f')
	if not options.summary_only:
		# Write the header descriptions
		dout.describeHeader(data_outs)
		# Write the header fields
		dout.writeHeader(data_outs)
	#format = dout.getFormat(named=True)
	n_written = 0

	sumout = util.DelimitedOutput()
	n_summary_written = 0
	if summary_output:
		sumout.addHeader('experiment.id','experiment identifier','s')
		sumout.addHeader('spin.rcf.g','rotor centrifugal force in g','e')
		sumout.addHeader('spin.time.sec','time of centrifuge spin in seconds','f')
		sumout.addHeader('crosslink.mol.prob','probability of crosslink per molecules','f')
		sumout.addHeader('crosslink.nt.prob','probability of crosslink per nucleotide','e')
		sumout.addHeader('translation.prob','fraction of maximal translation rate','e')
		sumout.addHeader('mol.translation.prob','probability this molecule is translated at maximal translation rate','e')
		sumout.addHeader('protection.prob','probability this molecule is protected from crosslinking','e')
		sumout.addHeader('id','molecule ID','s')
		sumout.addHeader('length','length in nucleotides','d')
		sumout.addHeader('mass','mass in g/mol','f')
		sumout.addHeader('radius','radius in nanometers','f')
		sumout.addHeader('psup','average proportion in the supernatant','f')
		sumout.addHeader('pfree','average proportion free (not in any cluster)','f')
		sumout.addHeader('psupclust','average proportion in the supernatant and in a cluster','f')
		sumout.addHeader('count','average number of molecules','f')
		sumout.addHeader('reps','number of replicates','d')
		sumout.addHeader('avg.cluster.count','average number of clusters this molecule is in','f')
		sumout.addHeader('avg.cluster.size','average number of molecules in average cluster','f')
		sumout.addHeader('avg.cluster.length','average cluster length in nucleotides','f')
		sumout.addHeader('avg.cluster.mass','average cluster mass in g/mol','e')
		sumout.addHeader('avg.cluster.radius','average cluster radius in nanometers','f')

		sumout.describeHeader(summary_outs)
		sumout.writeHeader(summary_outs)
		# If summary file requested, store data
		# Summary averages over molecules, then replicates, within experiments
		# Stores pSup
		# Key by (molecule_id, replicate_id)
		result_dict = {}

	for experiment_id in experiment_order:
		exp_params = experiment_settings_dict[experiment_id]

		# Set random seed
		np.random.seed(exp_params.random_seed)


		if summary_output:
			result_dict[experiment_id] = [] # refresh
		info_outs.write("\n# Experiment ID = {:s}".format(exp_params.experiment_id))
		info_outs.write("\n# Crosslink method = {:s}".format(exp_params.crosslink_method))
		info_outs.write("\n# Crosslink probability per molecule (e.g. 5' end) = {:1.2f}".format(exp_params.crosslink_mol_prob))
		info_outs.write("\n# Crosslink probability per nucleotide = {:1.2e}".format(exp_params.crosslink_nt_prob))
		try:
			trans_prob = exp_params.translation_prob
			info_outs.write("\n# Global proportion of maximal translation initiation rate = {:1.2e}\n".format(trans_prob))
		except KeyError:
			info_outs.write("\n")

		# Override options
		# DAD: implement
		#if not options.number_of_molecules is None:
		

		for rep in range(exp_params.replicates):
			info_outs.write("{:d} ".format(rep+1))
			info_outs.flush()

			if summary_output:
				rep_result_dict = {}

			# Generate molecules
			if not exp_params.template_filename is None:
				# Generate based on template transcriptome
				molecules, moldict, molecule_ids = sedsim.generateMoleculesFromTemplates(exp_params)
			else:
				# Generate de novo
				molecules, moldict, molecule_ids = sedsim.generateMolecules(exp_params.number_of_molecules)

			# Debugging
			#print(rep, exp_params.replicates)
			#mids = set([m.label for m in molecules])
			#info_outs.write("# Loaded {:d} molecules of {:d} types\n".format(len(molecules), len(mids)))

			# Generate replicates	
			clust_collection = sedsim.ClusterCollection()
			#xlink = sedsim.NucleotideCrosslinker(crosslink_avg)
			xlink = None
			if exp_params.crosslink_method == "DualCrosslinker":
				xlink = sedsim.DualCrosslinker(exp_params.crosslink_mol_prob, \
					exp_params.crosslink_nt_prob)
			if exp_params.crosslink_method == "DualVariableCrosslinker":
				xlink = sedsim.DualVariableCrosslinker(
					exp_params.crosslink_mol_prob, \
					exp_params.crosslink_nt_prob, \
					exp_params.translation_prob)

			# Add before crosslinking, giving us direct access to the input molecules
			for mol in molecules:
				clust_collection.add(mol)
			# Crosslink?
			if not xlink is None:
				clust_collection = xlink.crosslink(clust_collection)
			#print(len(clust_collection))
			tube = sedsim.SedimentationVessel(exp_params.sample_height_mm)
			#tube = sedsim.DeterministicSedimentationVessel(exp_params.sample_height_mm)
			tube.viscosity = exp_params.viscosity_mpa_sec
			tube.load(clust_collection)

			# Spin
			tube.spin(exp_params.rcf_g, exp_params.rotor_radius_cm, exp_params.spin_time_sec)

			mols = [m for m in tube.molecules]
			#print(mols)
			# Sort by size of cluster
			mols.sort(key=lambda x: x.container.count, reverse=True)
			num_clusters = 0
			for mol in mols:
				# A dictionary of results, one result per addHeader call above
				cluster = mol.container
				if cluster.label is None:
					num_clusters += 1
					cluster.label = "{:d}".format(num_clusters)
				if not options.summary_only:
					result = dout.createResult(default=None)
					result['replicate'] = rep+1
					result['experiment.id'] = exp_params.experiment_id
					result['spin.rcf.g'] = exp_params.rcf_g
					result['spin.time.sec'] = exp_params.spin_time_sec
					result['crosslink.mol.prob'] = exp_params.crosslink_mol_prob
					result['crosslink.nt.prob'] = exp_params.crosslink_nt_prob
					result['translation.prob'] = exp_params.translation.prob
					result['mol.translation.prob'] = mol.translation_prob
					result['mol.protection.prob'] = mol.protection_prob
					result['id'] = mol.label
					result['length'] = mol.length
					result['mass'] = int(mol.mass * 6.022e23)
					result['radius'] = mol.radius*1e9
					result['psup'] = cluster.psup
					result['cluster.id'] = cluster.label
					result['cluster.size'] = cluster.count
					result['cluster.length'] = cluster.length
					result['cluster.mass'] = int(cluster.mass * 6.022e23)
					result['cluster.radius'] = cluster.radius*1e9
					result['cluster.initial.position'] = cluster.initial_position_mm
					result['cluster.position'] = cluster.position_mm
					result['num.pelleted'] = 0
					if cluster.pelleted:
						result['num.pelleted'] = 1
					result['num.total'] = 1

				if summary_output:
					# Accumulate these results in the replicate storage dictionary
					if mol.label in rep_result_dict:
						rep_result_dict[mol.label].append(mol)
					else:
						rep_result_dict[mol.label] = [mol]
				#print(result)
				if not options.summary_only:
					# Parse the values, convert Nones to NA, etc.
					line = dout.formatLine(result)
					# A more manual approach:
					# line = format.format(column1="the answer is", column2=42)
					data_outs.write(line)
					n_written += 1

			# Now summarize molecules by molecule ID
			if summary_output:
				for mol_id in rep_result_dict.keys(): # mol_id = molecule ID
					molsum = sedsim.MoleculeSummary()
					mols = rep_result_dict[mol_id]
					#mols = [m for (m,c) in molresults]
					#clusters = [c for (m,c) in molresults]
					sentinel_mol = mols[0]
					molsum.id = mol_id
					molsum.length = sentinel_mol.length
					molsum.radius = sentinel_mol.radius
					molsum.mass = sentinel_mol.mass
					molsum.translation_prob = sentinel_mol.translation_prob
					molsum.count = len(mols)
					# psup, proportion in the supernatant
					molsum.psup = np.mean([m.psup for m in mols])
					# pfree, proportion free (not in clusters)
					molsum.pfree = np.mean([m.isFree() for m in mols])
					# psupclust, proportion in supernatant despite being in clusters 
					clustered_and_pelleted = [m.psup for m in mols if not m.isFree() and m.container.pelleted]
					if len(clustered_and_pelleted)>0:
						molsum.psupclust = np.mean(clustered_and_pelleted)
					else:
						molsum.psupclust = 0.0
					cc = [m.container.count for m in mols if m.container != m]
					if len(cc)>0:
						molsum.cluster_count = np.mean(cc)
					else:
						molsum.cluster_count = 1 # just itself
					molsum.cluster_length = np.mean([m.container.length for m in mols])
					molsum.cluster_size = np.mean([m.container.count for m in mols])
					molsum.cluster_mass = np.mean([m.container.mass for m in mols])
					molsum.cluster_radius = np.mean([m.container.radius for m in mols])
					# Replace the list of molecules with a single summary object
					rep_result_dict[mol_id] = molsum
				result_dict[experiment_id].append(rep_result_dict)


		# Done with replicates. Now summarize replicates.
		if summary_output:
			# Get list of replicates, each a dictionary of summaries keyed by molecule ID
			exp_results = result_dict[experiment_id]
			for molid in molecule_ids:
				try:
					mols = [res[molid] for res in exp_results if molid in res]
					if len(mols) > 0:
						result = sumout.createResult(default=None)
						result['experiment.id'] = experiment_id
						result['spin.rcf.g'] = options.rcf_g
						result['spin.time.sec'] = options.spin_time_sec
						result['crosslink.mol.prob'] = exp_params.crosslink_mol_prob
						result['crosslink.nt.prob'] = exp_params.crosslink_nt_prob
						sentinel_mol = mols[0]
						assert sentinel_mol.id == molid
						result['translation.prob'] = exp_params.translation_prob
						result['mol.translation.prob'] = sentinel_mol.translation_prob
						result['id'] = sentinel_mol.id
						result['length'] = sentinel_mol.length
						result['mass'] = int(sentinel_mol.mass * 6.022e23)
						result['radius'] = sentinel_mol.radius * 1e9
						result['psup'] = np.mean([m.psup for m in mols])
						result['pfree'] = np.mean([m.pfree for m in mols])
						result['psupclust'] = np.mean([m.psupclust for m in mols])
						result['count'] = np.mean([m.count for m in mols])
						result['reps'] = len(mols)
						result['avg.cluster.count'] = np.mean([m.cluster_count for m in mols])
						result['avg.cluster.size'] = np.mean([m.cluster_size for m in mols])
						result['avg.cluster.length'] = np.mean([m.cluster_length for m in mols])
						result['avg.cluster.mass'] = int(np.mean([m.cluster_mass for m in mols]) * 6.022e23)
						result['avg.cluster.radius'] = np.mean([m.cluster_radius for m in mols]) * 1e9
						line = sumout.formatLine(result)
						summary_outs.write(line)
						n_summary_written += 1
				except KeyError as ke:
					print(ke)
					# Possible for a molecule type to be not present
					# given stochasticity
					pass


	# Write out stopping time
	info_outs.write("\n# Run finished {}\n".format(util.timestamp()))
	if not options.summary_only:
		data_outs.write("\n# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None and not options.summary_only:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()
	# Shut down output
	if summary_output:
		info_outs.write("# Wrote {} lines to {}\n".format(n_summary_written, options.summary_out_fname))
		summary_outf.close()

