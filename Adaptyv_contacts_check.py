
import Bio.PDB
import numpy as np
import sys, os, glob


usage = '\n\nUsage: python Adaptyv_contacts_check.py folder_of_pdbs contact_cutoff_distance output_chimera_attribute_file\n\n'

if len(sys.argv) != 4:
	print(usage)
	sys.exit()


folder_of_pdbs = sys.argv[1]
contact_cutoff_distance = float(sys.argv[2])
output_chimera_attribute_file = sys.argv[3]


parser = Bio.PDB.PDBParser(QUIET = True)

#Credit for contact map: https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/

def calc_residue_dist(residue_one, residue_two) :
	"""Returns the C-alpha distance between two residues"""
	diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
	return np.sqrt(np.sum(diff_vector * diff_vector))


def calc_dist_matrix(chain_one, chain_two) :
	"""Returns a matrix of C-alpha distances between two chains"""
	answer = np.zeros((len(chain_one), len(chain_two)), float)
	for row, residue_one in enumerate(chain_one) :
		for col, residue_two in enumerate(chain_two) :
			answer[row, col] = calc_residue_dist(residue_one, residue_two)
	return answer


def get_A_B_contact_map(path_to_pdb, cutoff_distance):
	structure = parser.get_structure('scratch', path_to_pdb)
	model = structure[0]
	dist_matrix = calc_dist_matrix(model['A'], model['B'])
	contact_map = dist_matrix < cutoff_distance
	return contact_map


def get_B_res_in_contact_with_A(pdbfile):
	ctm = get_A_B_contact_map(pdbfile, contact_cutoff_distance)
	ctmT = np.transpose(ctm) #Now rows = res in chain B
	B_contacts_list = [True in row for row in ctmT] #True/False for all res in EGFR. True if it's in contact with chain A
	return B_contacts_list


full_contacts_matrix = [get_B_res_in_contact_with_A(pdbpath) for pdbpath in glob.glob(os.path.join(folder_of_pdbs, '*.pdb'))]

num_pdbs = len(full_contacts_matrix)

full_contacts_matrix = np.array(full_contacts_matrix)

full_contacts_sum = full_contacts_matrix.sum(axis = 0)

full_contacts_sum_norm = full_contacts_sum / num_pdbs


#Next write file for mapping onto structure in Chimera:

header_lines = ['#\n', '#  Binder contact frequency to map onto EGFR\n', '#\n', '#  From Adaptyv Bio Protein Design Competition\n', '#\n', '#  Use this file to assign the attribute in Chimera with the\n', '#  Define Attribute tool or the command defattr.\n', '#\n', 'attribute: contactfreq\n', 'match mode: 1-to-1\n', 'recipient: residues\n']

data_lines = [f'	:{i + 1}	{full_contacts_sum_norm[i]}\n' for i in range(len(full_contacts_sum_norm))]


with open(output_chimera_attribute_file, 'w') as outfile:
	for line in header_lines:
		outfile.write(line)
	for line in data_lines:
		outfile.write(line)


#testing:
print(len(full_contacts_matrix))
print(len(full_contacts_matrix[0]))
print(list(full_contacts_matrix[0]))


