#!/usr/bin/python

import sys
import math
import re
import os
import datetime
from itertools import chain
import argparse

def _read_args(): # reading input arguments

	parser = argparse.ArgumentParser(description='Adds new residue info to the given forcefield directory.')
	parser.add_argument('-m', '--mol2', nargs=1, metavar='mol2_FILE', help='original .mol2 file')
	parser.add_argument('-f', '--ff', nargs=1, metavar='ForceField_DIR', help='forcefield directory to edit. It should be in the same directory as this python file.')
	parser.add_argument('-n', '--name', nargs='?', metavar='residue_name', help='name of the new residue. If not supplied, the name of .mol2 file will be used.')
	group = parser.add_mutually_exclusive_group()
	group.add_argument('-a', '--acpype', metavar='acpype_DIR', help='name of the residue.acpype directory. It should be located in the same directory.')
	group.add_argument('-c', '--charmmgui', metavar='charmmgui_DIR', help='name of the charmm-gui-**** directory. It should be located in the same directory.')
	args = parser.parse_args()

	amber_flag = False
	charmm_flag = False
	forcefield = args.ff[0]

	if args.acpype:
		amber_flag = True
		acpype_dirn = str(args.acpype)
		acpype_resn = acpype_dirn.split('.', 1)[0]
	if args.charmmgui:
		charmm_flag = True
		charmm_dirn = str(args.charmmgui)
		ff = forcefield.split('-', 1)[0]
	
	if amber_flag:
		if not forcefield.startswith('amber'):
			sys.exit('You entered the acpype result, but the forcefield you are modifying is not AMBER!')
		else:
			if os.path.isfile(acpype_dirn + '/' + acpype_resn + '_GMX.itp'):
				itpFile = open(acpype_dirn + '/' + acpype_resn + '_GMX.itp', 'r')
				prmFile = None
			else:
				sys.exit("No " + acpype_resn + "_GMX.itp file in the .acpype folder!")
	if charmm_flag:
		if not forcefield.startswith('charmm'):
			sys.exit('You entered the charmm-gui result, but the forcefield you are modifying is not CHARMM!')
		else:
			if os.path.isfile(charmm_dirn + '/gromacs/LIG.itp'):
				itpFile = open(charmm_dirn + '/gromacs/LIG.itp', 'r')
			else:
				sys.exit("No LIG.itp file in the "+charmm_dirn+"/gromacs folder!")

			if os.path.isfile(charmm_dirn + '/gromacs/'+ff+'.itp'):
				prmFile = open(charmm_dirn + '/gromacs/'+ff+'.itp', 'r')
			else:
				sys.exit("No "+ff+".itp file in the "+charmm_dirn+"/gromacs folder!")

	mol2 = args.mol2[0]
	
	if os.path.isfile(mol2):
		mol2File = open(mol2, 'r')
	else:
		sys.exit("No such file as " + mol2 +"!")

	if args.name:
		resn = str(args.name)
	else:
		resn = mol2.split('.', 1)[0]

	resn = resn.upper()

	if len(resn) > 4:
		sys.exit("Residue name is too long!")
	
	if not os.path.isdir(forcefield):
		sys.exit("No such directory as "+ forcefield + "!")	

	return resn, mol2File, itpFile, prmFile, forcefield, amber_flag


def _cap_remover(mol2File, itpFile, prmFile, amber_flag): # recognize ACE and NME groups and remove them

	###############################################################
	# Reading the heavy atom names and bond info from .mol2 file  #
	###############################################################

	lines = mol2File.readlines()
	mol2File.close()
	
	atoms = []
	bonds = []
	flag = 0

	for line in lines:
		line = line.strip()
		if len(line) >0:
			segments = line.split()
			if segments[0] == '@<TRIPOS>MOLECULE':
				flag = 1
				continue
			elif segments[0] == '@<TRIPOS>ATOM':
				flag = 2
				continue
			elif segments[0] == '@<TRIPOS>BOND':
				flag = 3
				continue
			elif segments[0] == '@<TRIPOS>SUBSTRUCTURE':
				flag = 4

			if flag == 1:
				if len(segments) > 1:
					n_atoms = int(segments[0])
					n_bonds = int(segments[1])
			elif flag == 2:
				atoms.append({'index':int(segments[0]), 'name':segments[1], 'type':segments[5], 'resi':int(segments[6])})
			elif flag == 3:
				bonds.append([int(segments[1]), int(segments[2])])

	if len(atoms) != n_atoms:
		sys.exit("Number of atoms not matching in the .mol2 file!")
	if len(bonds) != n_bonds:
		sys.exit("Number of bonds not matching in the .mol2 file!")


	residues = set()
	
	for atom in atoms:
		residues.add(atom['resi'])

	if len(residues) != 3:
		sys.exit("More or less than 3 residues found! Only ACE+UNK+NME residues are allowed")
	else:
		print('mol2 file successfully read.')

	#################################
	# Finding ACE and NME cap atoms #
	#################################

	residue_index = list(residues)
	residue_index.sort()

	ACE_index = residue_index[0]
	UNK_index = residue_index[1]
	NME_index = residue_index[2]
	
	NME_Hs = []
	Hs = []
	Cs = []

	for atom in atoms:
		if atom['resi'] == ACE_index:
			if atom['type'] == 'C.3':
				ACE_methyl_index = atom['index']
			elif atom['type'] == 'C.2':
				ACE_C_index = atom['index']
			elif atom['type'] == 'O.2':
				ACE_O_index = atom['index']

		elif atom['resi'] == NME_index:
			if atom['type'] == 'N.am':
				NME_N_index = atom['index']
			elif atom['type'] == 'C.3':
				NME_methyl_index = atom['index']
			elif atom['type'] == 'H':
				NME_Hs.append(atom['index'])

		else:
			if atom['type'] == 'H':
				Hs.append(atom['index'])
			elif atom['type'].startswith('C'):
				Cs.append(atom['index'])
	
	for bond in bonds:
		if bond[0] == ACE_C_index:
			if bond[1] != ACE_methyl_index and bond[1] != ACE_O_index:
				N_index = bond[1]
		elif bond[1] == ACE_C_index:
			if bond[0] != ACE_methyl_index and bond[0] != ACE_O_index:
				N_index = bond[0]

		elif bond[0] == NME_N_index:
			if bond[1] != NME_methyl_index and bond[1] not in NME_Hs:
				C_index = bond[1]
		elif bond[1] == NME_N_index:
			if bond[0] != NME_methyl_index and bond[0] not in NME_Hs:
				C_index = bond[0]


	for bond in bonds:
		if bond[0] == NME_N_index:
			if bond[1] != NME_methyl_index and bond[1] != C_index:
				NME_H_index = bond[1]
		elif bond[1] == NME_N_index:
			if bond[0] != NME_methyl_index and bond[0] != C_index:
				NME_H_index = bond[0]
		elif bond[0] == C_index:
			if bond[1] in Cs:
				CA_index = bond[1]
			elif bond[1] != NME_N_index:
				O_index = bond[1]
		elif bond[1] == C_index:
			if bond[0] in Cs:
				CA_index = bond[0]
			elif bond[0] != NME_N_index:
				O_index = bond[0]

	
	proline_flag = False
	H_index = 0	

	for bond in bonds:
		if bond[0] == N_index:
			if bond[1] != ACE_C_index and bond[1] != CA_index:
				if bond[1] not in Hs:
					proline_flag = True
				else:
					H_index = bond[1]
		elif bond[1] == N_index:
			if bond[0] != ACE_C_index and bond[0] != CA_index:
				if bond[0] not in Hs:
					proline_flag = True
				else:
					H_index = bond[0]

	###########################
	# Renaming hydrogen atoms #
	###########################
	
	roots = {}
	
	for bond in bonds:
		if bond[0] in Hs:
			if bond[1] in roots:
				roots[bond[1]].append(bond[0])
			else:
				roots[bond[1]] = [bond[0]]
		elif bond[1] in Hs:
			if bond[0] in roots:
				roots[bond[0]].append(bond[1])
			else:
				roots[bond[0]] = [bond[1]]


	for atom in atoms:
		if atom['index'] in roots:
			mult = len(roots[atom['index']])
			order = 1
			for atom2 in atoms:
				if atom2['index'] in roots[atom['index']]:
					if mult > 1:
						name = 'H'+atom['name'][1:]+str(order)
					elif len(atom['name']) >1:
						name = 'H'+atom['name'][1:]
					elif atom['name'] == 'N':
						if amber_flag:
							name = 'H'
						else:
							name = 'HN'
					atom2['name'] = name
					order += 1

	#################################################
	# Reading atomtype info from input .itp file... #
	#################################################

	lines = itpFile.readlines()

	counter = -1
	flag = 0

	pairs = []
	angles = []
	propers = []
	impropers = []

	for line in lines:
		line = line.strip()
		if len(line) > 0:
			segments = line.split()
			if segments[0] == '[':
				if segments[1] == 'atoms':
					counter = 0
					continue
				elif segments[1] == 'pairs':
					flag = 1
					continue
				elif segments[1] == 'angles':
					flag = 2
					continue
				elif segments[1] == 'dihedrals':
					if flag < 3:
						flag = 3
					else:
						flag = 4
					continue

			if flag == 0 and counter >= 0 and counter < n_atoms and segments[0] != ';':
				atoms[counter]['type'] = segments[1]
				atoms[counter]['charge'] = float(segments[6])
				atoms[counter]['mass'] = float(segments[7])
				counter += 1

			if flag == 1 and segments[0] != ';':
				pairs.append([int(segments[0]), int(segments[1]), int(segments[2])])

			if flag == 2 and segments[0] != ';':
				angles.append([int(segments[0]), int(segments[1]), int(segments[2])])

			if flag == 3 and segments[0] != ';':
				propers.append([int(segments[0]), int(segments[1]), int(segments[2]), int(segments[3])])

			if flag == 4 and segments[0] != ';':
				impropers.append([int(segments[0]), int(segments[1]), int(segments[2]), int(segments[3])])


	##########################################
	# Removing charges from ACE and NME caps #
	##########################################


	cap_group = []
	cap_charge = 0.0
	total_charge = 0.0

	for atom in atoms:
		total_charge += atom['charge']
		if atom['resi'] == ACE_index or atom['resi'] == NME_index:
			cap_group.append(atom['index'])
			cap_charge += atom['charge']

	new_total_charge = float(round(total_charge))

	print(f"The given molecule has a total charge of {new_total_charge:.6f}")

	backbone = [N_index, H_index, C_index, O_index]

	newAtoms = []

	for atom in atoms:
		if atom['resi'] == UNK_index:
			atomType = atom['type']
			charge = atom['charge']
			mass = atom['mass']
			
			if atom['index'] == N_index:
				if not proline_flag:
					if amber_flag:
						atomType = 'N' 
					else:
						atomType = 'NH1'
						charge = -0.4700
				else:
					if amber_flag:
						atomType = 'N'
					else:
						atomType = 'N'
						charge = -0.2900
			elif atom['index'] == H_index:
				if not proline_flag:
					if amber_flag:
						atomType = 'H'
					else:
						atomType = 'H' 
						charge = 0.3100
			elif atom['index'] == C_index:
				if amber_flag:
					atomType = 'C'
				else:
					atomType = 'C' 
					charge = 0.5100
			elif atom['index'] == O_index:
				if amber_flag:
					atomType = 'O'
				else:
					atomType = 'O' 
					charge = -0.5100
			elif atom['index'] == CA_index:
				if amber_flag:
					atomType = 'CX'
				else:
					atomType = 'CT1'
				
			newAtoms.append({'index':atom['index'], 'name':atom['name'], 'type':atomType, 'charge':charge, 'mass':mass})

	total_charge = 0.0	

	for atom in newAtoms:
		total_charge += atom['charge']

	charge_change = new_total_charge - total_charge

	print(f"Total charge should be changed by {charge_change:.6f}.")

	if charge_change > 0:
		sign = 1
	else:
		sign = -1

	steps = int(abs(round(charge_change, 6)*1000000))

	print(f"Total {steps} steps will be run to compensate the cap charge.")

	i = 0

	while i < steps:
		resi = i%len(newAtoms)
		if amber_flag or newAtoms[resi]['index'] not in backbone:
			newAtoms[resi]['charge'] += 0.000001*sign
		else:
			steps += 1
		i += 1

	charge_sum = 0.0
	for atom in newAtoms:
		atom['charge'] = round(atom['charge'],6)
		charge_sum += atom['charge']

	if round(charge_sum, 6) != new_total_charge:
		print(f"Current charge: {charge_sum}")
		sys.exit("Error in charge compensation!")

	else:
		print("Charge redistribution complete.")


	####################################################################################
	# Relabeling bonds/pairs/angles/propers/impropers info from indices into atomtypes #
	####################################################################################
	
	index2type = {}

	for atom in newAtoms:
		index2type[atom['index']] = atom['type']

	for atom in atoms:
		if atom['resi'] == ACE_index:
			if atom['index'] == ACE_C_index:
				index2type[atom['index']] = 'C'
			elif atom['index'] == ACE_O_index:
				index2type[atom['index']] = 'O'
			elif atom['index'] == ACE_methyl_index:
				if amber_flag:
					index2type[atom['index']] = 'CX'
				else:
					index2type[atom['index']] = 'CT1'
			else:
				if amber_flag:
					index2type[atom['index']] = 'H1'
				else:
					index2type[atom['index']] = 'HB1'

		elif atom['resi'] == NME_index:
			if atom['index'] == NME_N_index:
				if amber_flag:
					index2type[atom['index']] = 'N'
				else:
					index2type[atom['index']] = 'NH1'
			elif atom['index'] == NME_H_index:
				index2type[atom['index']] = 'H'
			elif atom['index'] == NME_methyl_index:
				if amber_flag:
					index2type[atom['index']] = 'CX'
				else:
					index2type[atom['index']] = 'CT1'
			else:
				if amber_flag:
					index2type[atom['index']] = 'H1'
				else:
					index2type[atom['index']] = 'HB1'


	############################
	# Removing cap group infos #
	############################

	newBonds = []
	newPairs = []
	newAngles = []
	newPropers = []
	newImpropers = []

	for bond in bonds:
		if not (bond[0] in cap_group or bond[1] in cap_group):
			newBonds.append(bond)

	for pair in pairs:
		if not (pair[0] in cap_group and pair[1] in cap_group):
			newPairs.append(pair)
	
	for angle in angles:
		if not ((angle[0] in cap_group and angle[1] in cap_group) or (angle[1] in cap_group and angle[2] in cap_group)):
			newAngles.append(angle)

	for proper in propers:
		if not ((proper[0] in cap_group and proper[1] in cap_group and proper[2] in cap_group) or (proper[1] in cap_group and proper[2] in cap_group and proper[3] in cap_group)):
			newPropers.append(proper)

	for improper in impropers:
		if (amber_flag and not improper[2] in cap_group) or (not amber_flag and not improper[0] in cap_group):
			newImpropers.append(improper)

	Impropers = [] # Improper info for .rtp file

	for improper in impropers:
		if not improper[0] in cap_group and not improper[1] in cap_group and not improper[2] in cap_group and not improper[3] in cap_group:
			Impropers.append(improper)



	####################################################
	# Reading parameter values from .itp or .prm files #
	####################################################

	atomList = []
	bondList = []
	pairList = []
	angleList = []
	properList = []
	improperList = []

	if amber_flag:
		itpFile.seek(0, 0)
		lines = itpFile.readlines()
		itpFile.close()

		flag = 0
		
		for line in lines:
			line = line.strip()
			if len(line) > 0:
				segments = line.split()
				if segments[0] != ';' and not segments[0].startswith(';'):
					if segments[0] == '[':
						if segments[1] == 'atomtypes':
							flag = 1
							continue
						elif segments[1] == 'moleculetype':
							flag = 2
							continue
						elif segments[1] == 'atoms':
							flag = 3
							continue
						elif segments[1] == 'bonds':
							flag = 4
							continue
						elif segments[1] == 'pairs':
							flag = 5
							continue
						elif segments[1] == 'angles':
							flag = 6
							continue
						elif segments[1] == 'dihedrals':
							if flag == 6:
								flag = 7
								continue
							else:
								flag = 8
								continue
					if flag == 1: # atomtypes for AMBER forcefield
						name = segments[0]

						if name.startswith(';'):
							continue

						elif name == 'cl':
							anum = 17
							mass = 35.453
						elif name == 'br':
							anum = 35
							mass = 79.904
						elif name == 'i':
							anum = 53
							mass = 126.904
						elif name.upper().startswith('C'):
							anum = 6
							mass = 12.01
						elif name.upper().startswith('H'):
							anum = 1
							mass = 1.008
						elif name.upper().startswith('N'):
							anum = 7
							mass = 14.01
						elif name.upper().startswith('O'):
							anum = 8
							mass = 16.00
						elif name.upper().startswith('F'):
							anum = 9
							mass = 19.00
						elif name.upper().startswith('P'):
							anum = 15
							mass = 30.97
						elif name.upper().startswith('S'):
							anum = 16
							mass = 32.06
						else:
							sys.exit("Unknown atomtype found in the GMX.itp file!")

						ptype = segments[4]
						sigma = segments[5]
						epsilon = segments[6]

						atomList.append({'name':name, 'anum':anum, 'mass':mass, 'ptype':ptype, 'sigma':sigma, 'epsilon':epsilon})
					elif flag == 2: # skipping moleculetype
						continue
					elif flag == 3: # skipping atoms
						continue

					elif flag == 4: # bondtypes for AMBER forcefield
						i = int(segments[0])
						j = int(segments[1])
						func = segments[2]
						b0 = float(segments[3])
						kb = float(segments[4])
								
						if [i,j] in newBonds:
							i = index2type[i]
							j = index2type[j]
							bondList.append({'i':i, 'j':j, 'func':func, 'b0':b0, 'kb':kb})
					elif flag == 5: # skipping pairtypes
						continue
					elif flag == 6: # angletypes for AMBER forcefield
						i = int(segments[0])
						j = int(segments[1])
						k = int(segments[2])
						func = segments[3]
						th0 = float(segments[4])
						cth = float(segments[5])

						if [i,j,k] in newAngles:
							i = index2type[i]
							j = index2type[j]
							k = index2type[k]
							angleList.append({'i':i, 'j':j, 'k':k, 'func':func, 'th0':th0, 'cth':cth})

					elif flag == 7: # proper dihedral parameters
						i = int(segments[0])
						j = int(segments[1])
						k = int(segments[2])
						l = int(segments[3])
						func = segments[4] 
						phase = float(segments[5]) 
						kd = float(segments[6])
						pn = segments[7]

						if [i,j,k,l] in newPropers:
							i = index2type[i]
							j = index2type[j]
							k = index2type[k]
							l = index2type[l]
							properList.append({'i':i, 'j':j, 'k':k, 'l':l, 'func':func, 'phase':phase, 'kd':kd, 'pn':pn})
					elif flag == 8: # improper dihedral parameters
						i = int(segments[0])
						j = int(segments[1])
						k = int(segments[2])
						l = int(segments[3])
						func = segments[4] 
						phase = float(segments[5]) 
						kd = float(segments[6])
						pn = segments[7]

						if [i,j,k,l] in newImpropers:
							i = index2type[i]
							j = index2type[j]
							k = index2type[k]
							l = index2type[l]
							improperList.append({'i':i, 'j':j, 'k':k, 'l':l, 'func':func, 'phase':phase, 'kd':kd, 'pn':pn})
		

	else: # for CHARMM forcefield
		lines = prmFile.readlines()
		prmFile.close()

		flag = 0
		
		for line in lines:
			line = line.strip()
			if len(line) > 0:
				segments = line.split()
				if segments[0] != ';':
					if segments[0] == '[':
						if segments[1] == 'atomtypes':
							flag = 1
							continue
						elif segments[1] == 'bondtypes':
							flag = 2
							continue
						elif segments[1] == 'pairtypes':
							flag = 3
							continue
						elif segments[1] == 'angletypes':
							flag = 4
							continue
						elif segments[1] == 'dihedraltypes':
							if flag == 4:
								flag = 5
								continue
							else:
								flag = 6
								continue

					if flag == 1: 
						atomList.append({'name':segments[0], 'anum':segments[1], 'mass':float(segments[2]), 'charge':float(segments[3]), 'ptype':segments[4], 'sigma':float(segments[5]), 'epsilon':float(segments[6])})
					elif flag == 2:
						bondList.append({'i':segments[0], 'j':segments[1], 'func':segments[2], 'b0':float(segments[3]), 'kb':float(segments[4])})
					elif flag == 3:
						pairList.append({'i':segments[0], 'j':segments[1], 'func':segments[2], 'sig':float(segments[3]), 'eps':float(segments[4])})
					elif flag == 4:
						angleList.append({'i':segments[0], 'j':segments[1], 'k':segments[2], 'func':segments[3], 'th0':float(segments[4]), 'cth':float(segments[5]), 'S0':float(segments[6]), 'Kub':float(segments[7])})
					elif flag == 5:
						properList.append({'i':segments[0], 'j':segments[1], 'k':segments[2], 'l':segments[3], 'func':segments[4], 'phase':float(segments[5]), 'kd':float(segments[6]), 'pn':segments[7]})
					elif flag == 6:
						improperList.append({'i':segments[0], 'j':segments[1], 'k':segments[2], 'l':segments[3], 'func':segments[4], 'phase':float(segments[5]), 'kd':float(segments[6])})


	return atoms, newAtoms, bonds, newBonds, newPairs, newAngles, newPropers, newImpropers, Impropers, atomList, bondList, pairList, angleList, properList, improperList, roots



def _rtp_updator(resn, atoms, bonds, impropers, forcefield, amber_flag): # appends new residue info to the aminoacids.rtp file
	os.chdir(forcefield)
	print('Appending to ../'+forcefield +'/aminoacids.rtp file...')
	
	rtpFile = open('aminoacids.rtp', 'a')
	
	now = datetime.datetime.now() 
	nowDate = now.strftime('%Y-%m-%d')

	
	"""
	In the aminoacids.rtp file, the bonds/improper info must be written in atom names!
	"""	

	rtpFile.write('\n\n; '+resn+' added by '+nowDate+'\n\n')
	rtpFile.write('\n[ '+resn+' ]'+'\n')
	rtpFile.write('  [ atoms ]\n')

	groupCounter = 0

	if amber_flag:
		for atom in atoms:
			groupCounter += 1
			rtpFile.write(f"{atom['name']:>6s}{atom['type']:>6s}{atom['charge']:>19.6f}{groupCounter:>6d}\n")
	else:
		for atom in atoms:
			groupCounter += 1
			rtpFile.write(f"{atom['name']:>9s}{atom['type']:>9s}{atom['charge']:>8.6f}{groupCounter:>3d}\n")

	rtpFile.write('  [ bonds ]\n')

	for bond in bonds:
		for atom in atoms:
			if atom['index'] == bond[0]:
				a1 = atom['name']
			if atom['index'] == bond[1]:
				a2 = atom['name']
		if amber_flag:
			rtpFile.write(f"{a1:>6s}{a2:>6s}\n")
		else:
			rtpFile.write(f"{a1:>9s}{a2:>6s}\n")

	if amber_flag:
		rtpFile.write("    -C     N\n")
	else:
		rtpFile.write("        C    +N\n")

	rtpFile.write('  [ impropers ]\n')

	if amber_flag:
		rtpFile.write('    -C    CA     N     H\n')
		rtpFile.write('    CA    +N     C     O\n')
	else:
		rtpFile.write('        N    -C    CA    HN\n')
		rtpFile.write('        C    CA    +N     O\n')
	

	for improper in impropers:
		for atom in atoms:
			if atom['index'] == improper[0]:
				a1 = atom['name']
			if atom['index'] == improper[1]:
				a2 = atom['name']
			if atom['index'] == improper[2]:
				a3 = atom['name']
			if atom['index'] == improper[3]:
				a4 = atom['name']
		if amber_flag:
			rtpFile.write(f"{a1:>6s}{a2:>6s}{a3:>6s}{a4:>6s}\n")
		else:
			rtpFile.write(f"{a1:>9s}{a2:>6s}{a3:>6s}{a4:>6s}\n")

	if not amber_flag:
		rtpFile.write('  [ cmap ]\n')
		rtpFile.write('       -C     N    CA     C    +N\n\n')

	rtpFile.close()	

	return 0


def _hydrogen_sorter(atoms, bonds, hydrogens, amber_flag): # group hydrogens and clarify root atoms

	#-----------------------------
	# AMBER hydrogen type library
	#-----------------------------

	group_1 = {'CA','CB','CC','CK','CM','CN','CQ','CR','CV','CW','C*','N','NA','NB','NC','ca','c','n','nh','n2','cc','cd','ce','cf','cu','cv','nb','nc','nd'}
	group_2 = {'OH','O2','o','oh','sh'}
	group_3 = {'N','N2','c','c2','n','nh'}
	group_4 = {'CX','CT','n4','c3'}
	group_5 = {'CX','CI','CJ','CT','N3','n4','c3','n3'}
	group_6 = {'CX','CT','c3','cx','cy'}


	#------------------------------
	# CHARMM hydrogen type library
	#------------------------------

	group1 = {'H','HL','HGR51','HGR52','HGR53','HGR61','HGR62','HGR63','HGR71','HGP1','HGP2','HN2','HN3','HP','HR1','HR2','HR3','HNP','HN3B'}
	group2 = {'HCP1','HCP1M','HGP2','HGP3','HN4','HN5','HOL','HS','HSIO','HPER','HX'}
	group3 = {'HGA4','HGA5','HGP4','HEL1','HEL2','HN1','HE1','HE2','HNE1','HNE2'}
	group4 = {'HB','HCA2A','HCA3A','HCA3','HCA3M','HGA3','HGAAM0','HGAAM1','HGAAM2','HGP5','HGPAM1','HGPAM2','HGPAM3','HAL3','HBL','HN6','HN8','HN9','HA3','HC','HAL','HSIA','HRA3','HRM1','HRM2'}
	group5 = {'HCA1A','HCA1','HCA1C2','HGA1','HAL1','HA','HA1','HB1'}
	group6 = {'HF1','HF2','HCA25A','HCA2','HCA2C2','HGA2','HGA6','HGA7','HAL2','HA','HA2','HB2'}
	group7 = {'HWT3','HWSPCE'}
	group8 = {'HWT4','HWT4EW'}
	group9 = {'HWT5'}


	n_groups = len(hydrogens)
	Hs = list(chain.from_iterable(list(hydrogens.values())))
	roots = list(hydrogens.keys())

	H_group = []
	
	for atom in atoms:
		if atom['index'] in roots:
			mult = len(hydrogens[atom['index']])
			neighbor = []
			for bond in bonds:
				if atom['index'] == bond[0]:
					if bond[1] not in Hs:
						neighbor.append(bond[1])
				elif atom['index'] == bond[1]:
					if bond[0] not in Hs:
						neighbor.append(bond[0])
			neighbor.sort()

			if amber_flag:
				if mult == 1:
					if atom['type'] in group_1:
						grp = 1
					elif atom['type'] in group_2:
						grp = 2
					elif atom['type'] in group_5:
						grp = 5
					else:
						sys.exit("Error!")
				elif mult == 2:
					if atom['type'] in group_3:
						grp = 3
					elif atom['type'] in group_4:
						grp = 4
					elif atom['type'] in group_6:
						grp = 6
				elif mult == 3:
					if atom['type'] in group_4:
						grp = 4

			else:
				H = hydrogens[atom['index']][0]
				for atom2 in atoms:
					if atom2['index'] == H:
						if atom2['type'] in group1:
							grp = 1
						elif atom2['type'] in group2:
							grp = 2
						elif atom2['type'] in group3:
							grp = 3
						elif atom2['type'] in group4:
							grp = 4
						elif atom2['type'] in group5:
							grp = 5
						elif atom2['type'] in group6:
							grp = 6
						elif atom2['type'] in group7:
							grp = 7
						elif atom2['type'] in group8:
							grp = 8
						elif atom2['type'] in group9:
							grp = 9
			
			atomA = ''
			atomB = ''
			atomC = ''			

			for i, index in enumerate(neighbor):
				for atom2 in atoms:
					if index == atom2['index']:
						if i == 0:
							atomA = atom2['name']
							atomAi = atom2['index']
						elif i == 1:
							atomB = atom2['name']
						elif i == 2:
							atomC = atom2['name']

			atomBi = -1


			if atomB == '':
				for bond in bonds:
					if bond[0] == atomAi and bond[1] not in Hs and bond[1] != atom['index']:
						atomBi = bond[1]
					elif bond[1] == atomAi and bond[0] not in Hs and bond[0] != atom['index']:
						atomBi = bond[0]

			for atom2 in atoms:
				if atom2['index'] == atomBi:
					atomB = atom2['name']
		
			H_group.append({'mult':mult, 'grp':grp, 'root':atom['name'], 'A':atomA, 'B':atomB, 'C':atomC})

	return H_group


def _hdb_updator(resn, H_group, forcefield, amber_flag): # appends new residue info to the aminoacids.hdb file

	print('Appending to ../'+forcefield +'/aminoacids.hdb file...')
	
	hdbFile = open('aminoacids.hdb', 'a')

	hdbFile.write(resn + '\t\t'+ str(len(H_group))+'\n')

	for H_grp in H_group:
		if H_grp['root'] == 'N':
			if amber_flag:
				hdbFile.write("1\t1\tH\tN\t-C\tCA\n")
			else:
				hdbFile.write("1\t1\tHN\tN\tCA\t-C\n")
		else:
			hdbFile.write(f"{H_grp['mult']}\t{H_grp['grp']}\tH{H_grp['root'][1:]}\t{H_grp['root']}\t{H_grp['A']}\t{H_grp['B']}\t{H_grp['C']}\n")

	hdbFile.close()
	
	return 0


def _list_updator(oldAtoms, atoms, bonds, pairs, angles, propers, impropers, bondList, pairList, angleList, properList, improperList, amber_flag):
	if amber_flag:
		return bondList, pairList, angleList, properList, improperList

	else:
		totalBondList = []
		totalPairList = []
		totalAngleList = []
		totalProperList = []
		totalImproperList = []

		for bond in bonds:
			for atom in oldAtoms:
				if atom['index'] == bond[0]:
					i = atom['index']
					ai = atom['type']
				elif atom['index'] == bond[1]:
					j = atom['index']
					aj = atom['type']
			totalBondList.append({'i':i, 'j':j, 'ai':ai, 'aj':aj, 'func':0, 'b0':0.0, 'kb':0.0})

		for pair in pairs:
			for atom in oldAtoms:
				if atom['index'] == pair[0]:
					i = atom['index']
					ai = atom['type']
				elif atom['index'] == pair[1]:
					j = atom['index']
					aj = atom['type']
			totalPairList.append({'i':i, 'j':j, 'ai':ai, 'aj':aj, 'func':0, 'sig':0.0, 'eps':0.0})

		
		for angle in angles:
			for atom in oldAtoms:
				if atom['index'] == angle[0]:
					i = atom['index']
					ai = atom['type']
				elif atom['index'] == angle[1]:
					j = atom['index']
					aj = atom['type']
				elif atom['index'] == angle[2]:
					k = atom['index']
					ak = atom['type']
			totalAngleList.append({'i':i, 'j':j, 'k':k, 'ai':ai, 'aj':aj, 'ak':ak, 'func':0, 'th0':0.0, 'cth':0.0, 'S0':0.0, 'Kub':0.0})

	
		for proper in propers:
			for atom in oldAtoms:
				if atom['index'] == proper[0]:
					i = atom['index']
					ai = atom['type']
				elif atom['index'] == proper[1]:
					j = atom['index']
					aj = atom['type']
				elif atom['index'] == proper[2]:
					k = atom['index']
					ak = atom['type']
				elif atom['index'] == proper[3]:
					l = atom['index']
					al = atom['type']
			totalProperList.append({'i':i, 'j':j, 'k':k, 'l':l, 'ai':ai, 'aj':aj, 'ak':ak, 'al':al, 'func':0, 'phase':0.0, 'kd':0.0, 'pn':0})

	
		for improper in impropers:
			for atom in oldAtoms:
				if atom['index'] == improper[0]:
					i = atom['index']
					ai = atom['type']
				elif atom['index'] == improper[1]:
					j = atom['index']
					aj = atom['type']
				elif atom['index'] == improper[2]:
					k = atom['index']
					ak = atom['type']
				elif atom['index'] == improper[3]:
					l = atom['index']
					al = atom['type']
			totalImproperList.append({'i':i, 'j':j, 'k':k, 'l':l, 'ai':ai, 'aj':aj, 'ak':ak, 'al':al, 'func':0, 'phase':0.0, 'kd':0.0})

		
		for bond in totalBondList:
			bonFile = open('ffbonded.itp', 'r')

			flag = 0

			lines = bonFile.readlines()
			bonFile.close()
			
			for line in lines:
				line = line.strip()
				if len(line) > 0:
					segments = line.split()
					if segments[0] != ';':
						if segments[0] == '[':
							if segments[1] == 'bondtypes':
								flag = 1
								continue
							elif segments[1] == 'angletypes':
								flag = 2
								continue
						if flag == 1:
							if bond['ai'] == segments[0] and bond['aj'] == segments[1]:
								bond['func'] = segments[2]
								bond['b0'] = float(segments[3])
								bond['kb'] = float(segments[4])		
								break
							if bond['ai'] == segments[1] and bond['aj'] == segments[0]:
								bond['func'] = segments[2]
								bond['b0'] = float(segments[3])
								bond['kb'] = float(segments[4])					
								break
						if flag > 1:
							break	

			if bond['func'] == 0:
				for bond2 in bondList:
					if bond['ai'] == bond2['i'] and bond['aj'] == bond2['j']:
						bond['func'] = bond2['func']
						bond['b0'] = bond2['b0']
						bond['kb'] = bond2['kb']
						break
					if bond['ai'] == bond2['j'] and bond['aj'] == bond2['i']:
						bond['func'] = bond2['func']
						bond['b0'] = bond2['b0']
						bond['kb'] = bond2['kb']			
						break

			if bond['func'] == 0:
				sys.exit("Error: bond info not found!")


		for pair in totalPairList:
			nonbonFile = open('ffnonbonded.itp', 'r')

			flag = 0

			lines = nonbonFile.readlines()
			nonbonFile.close()
			
			for line in lines:
				line = line.strip()
				if len(line) > 0:
					segments = line.split()
					if segments[0] != ';':
						if segments[0] == '[':
							if segments[1] == 'atomtypes':
								flag = 1
								continue
							elif segments[1] == 'pairtypes':
								flag = 2
								continue
						if flag == 2:
							if pair['ai'] == segments[0] and pair['aj'] == segments[1]:
								pair['func'] = segments[2]
								pair['sig'] = float(segments[3])
								pair['eps'] = float(segments[4])		
								break
							if pair['ai'] == segments[1] and pair['aj'] == segments[0]:
								pair['func'] = segments[2]
								pair['sig'] = float(segments[3])
								pair['eps'] = float(segments[4])		
								break

			if pair['func'] == 0:
				for pair2 in pairList:
					if pair['ai'] == pair2['i'] and pair['aj'] == pair2['j']:
						pair['func'] = pair2['func']
						pair['sig'] = pair2['sig']
						pair['eps'] = pair2['eps']
						break
					if pair['ai'] == pair2['j'] and pair['aj'] == pair2['i']:
						pair['func'] = pair2['func']
						pair['sig'] = pair2['sig']
						pair['eps'] = pair2['eps']			
						break

			if pair['func'] == 0:
				sys.exit("Error: pair info not found!")

		
		for angle in totalAngleList:
			bonFile = open('ffbonded.itp', 'r')

			flag = 0

			lines = bonFile.readlines()
			bonFile.close()
			
			for line in lines:
				line = line.strip()
				if len(line) > 0:
					segments = line.split()
					if segments[0] != ';':
						if segments[0] == '[':
							if segments[1] == 'angletypes':
								flag = 1
								continue
							elif segments[1] == 'dihedraltypes':
								flag = 2
								continue
						if flag == 1:
							if angle['ai'] == segments[0] and angle['aj'] == segments[1] and angle['ak'] == segments[2]:
								angle['func'] = segments[3]
								angle['th0'] = float(segments[4])
								angle['cth'] = float(segments[5])
								angle['S0'] = float(segments[6])
								angle['Kub'] = float(segments[7])		
								break
							if angle['ai'] == segments[2] and angle['aj'] == segments[1] and angle['ak'] == segments[0]:
								angle['func'] = segments[3]
								angle['th0'] = float(segments[4])
								angle['cth'] = float(segments[5])
								angle['S0'] = float(segments[6])
								angle['Kub'] = float(segments[7])		
								break
						if flag > 1:
							break	

			if angle['func'] == 0:
				for angle2 in angleList:
					if angle['ai'] == angle2['i'] and angle['aj'] == angle2['j'] and angle['ak'] == angle2['k']:
						angle['func'] = angle2['func']
						angle['th0'] = angle2['th0']
						angle['cth'] = angle2['cth']
						angle['S0'] = angle2['S0']
						angle['Kub'] = angle2['Kub']
						break
					if angle['ai'] == angle2['k'] and angle['aj'] == angle2['j'] and angle['ak'] == angle2['i']:
						angle['func'] = angle2['func']
						angle['th0'] = angle2['th0']
						angle['cth'] = angle2['cth']
						angle['S0'] = angle2['S0']
						angle['Kub'] = angle2['Kub']
						break

			if angle['func'] == 0:
				sys.exit("Error: Angle info not found!")


		for proper in totalProperList:
			bonFile = open('ffbonded.itp', 'r')

			flag = 0

			lines = bonFile.readlines()
			bonFile.close()
			
			for line in lines:
				line = line.strip()
				if len(line) > 0:
					segments = line.split()
					if segments[0] != ';':
						if segments[0] == '[':
							if segments[1] == 'dihedraltypes' and flag == 0:
								flag = 1
								continue
							elif segments[1] == 'dihedraltypes' and flag == 1:
								flag = 2
								continue
						if flag == 1:
							if proper['ai'] == segments[0] and proper['aj'] == segments[1] and proper['ak'] == segments[2] and proper['al'] == segments[3]:
								proper['func'] = segments[4]
								proper['phase'] = float(segments[5])
								proper['kd'] = float(segments[6])
								proper['pn'] = segments[7]		
								break
							if proper['ai'] == segments[3] and proper['aj'] == segments[2] and proper['ak'] == segments[1] and proper['al'] == segments[0]:
								proper['func'] = segments[4]
								proper['phase'] = float(segments[5])
								proper['kd'] = float(segments[6])
								proper['pn'] = segments[7]		
								break
						if flag > 1:
							break	

			if proper['func'] == 0:
				for proper2 in properList:
					if proper['ai'] == proper2['i'] and proper['aj'] == proper2['j'] and proper['ak'] == proper2['k'] and proper['al'] == proper2['l']:
						proper['func'] = proper2['func']
						proper['phase'] = proper2['phase']
						proper['kd'] = proper2['kd']
						proper['pn'] = proper2['pn']
						break
					if proper['ai'] == proper2['l'] and proper['aj'] == proper2['k'] and proper['ak'] == proper2['j'] and proper['al'] == proper2['i']:
						proper['func'] = proper2['func']
						proper['phase'] = proper2['phase']
						proper['kd'] = proper2['kd']
						proper['pn'] = proper2['pn']
						break

			if proper['func'] == 0:
				sys.exit("Error: proper dihedral info not found!")


		for improper in totalImproperList:
			bonFile = open('ffbonded.itp', 'r')

			flag = 0

			lines = bonFile.readlines()
			bonFile.close()
			
			for line in lines:
				line = line.strip()
				if len(line) > 0:
					segments = line.split()
					if segments[0] != ';':
						if segments[0] == '[':
							if segments[1] == 'dihedraltypes' and flag == 0:
								flag = 1
								continue
							elif segments[1] == 'dihedraltypes' and flag == 1:
								flag = 2
								continue
						if flag == 2:
							if improper['ai'] == segments[0] and improper['aj'] == segments[1] and improper['ak'] == segments[2] and improper['al'] == segments[3]:
								improper['func'] = segments[4]
								improper['phase'] = float(segments[5])
								improper['kd'] = float(segments[6])	
								break
							if improper['ai'] == segments[3] and improper['aj'] == segments[2] and improper['ak'] == segments[1] and improper['al'] == segments[0]:
								improper['func'] = segments[4]
								improper['phase'] = float(segments[5])
								improper['kd'] = float(segments[6])	
								break

			if improper['func'] == 0:
				for improper2 in improperList:
					if improper['ai'] == improper2['i'] and improper['aj'] == improper2['j'] and improper['ak'] == improper2['k'] and improper['al'] == improper2['l']:
						improper['func'] = improper2['func']
						improper['phase'] = improper2['phase']
						improper['kd'] = improper2['kd']
						break
					if improper['ai'] == improper2['l'] and improper['aj'] == improper2['k'] and improper['ak'] == improper2['j'] and improper['al'] == improper2['i']:
						improper['func'] = improper2['func']
						improper['phase'] = improper2['phase']
						improper['kd'] = improper2['kd']
						improper['pn'] = improper2['pn']
						break

			if improper['func'] == 0:
				sys.exit("Error: improper dihedral info not found!")


		newBondList = []
		newPairList = []
		newAngleList = []
		newProperList = []
		newImproperList = []

		for bond in totalBondList:
			for atom in atoms:
				if atom['index'] == bond['i']:
					i = atom['type']
				elif atom['index'] == bond['j']:
					j = atom['type']
			newBondList.append({'i':i,'j':j, 'func': bond['func'], 'b0': bond['b0'], 'kb': bond['kb']})

		for pair in totalPairList:
			for atom in atoms:
				if atom['index'] == pair['i']:
					i = atom['type']
				elif atom['index'] == pair['j']:
					j = atom['type']
			newPairList.append({'i':i, 'j':j, 'func':pair['func'], 'sig':pair['sig'], 'eps':pair['eps']})

		for angle in totalAngleList:
			for atom in atoms:
				if atom['index'] == angle['i']:
					i = atom['type']
				elif atom['index'] == angle['j']:
					j = atom['type']
				elif atom['index'] == angle['k']:
					k = atom['type']
			newAngleList.append({'i':i, 'j':j, 'k':k, 'func':angle['func'], 'th0':angle['th0'], 'cth':angle['cth'], 'S0':angle['S0'], 'Kub':angle['Kub']})

		for proper in totalProperList:
			for atom in atoms:
				if atom['index'] == proper['i']:
					i = atom['type']
				elif atom['index'] == proper['j']:
					j = atom['type']
				elif atom['index'] == proper['k']:
					k = atom['type']
				elif atom['index'] == proper['l']:
					l = atom['type']
			newProperList.append({'i':i, 'j':j, 'k':k, 'l':l, 'func':proper['func'], 'phase':proper['phase'], 'kd':proper['kd'], 'pn':proper['pn']})


		for improper in totalImroperList:
			for atom in atoms:
				if atom['index'] == improper['i']:
					i = atom['type']
				elif atom['index'] == improper['j']:
					j = atom['type']
				elif atom['index'] == improper['k']:
					k = atom['type']
				elif atom['index'] == improper['l']:
					l = atom['type']
			newImproperList.append({'i':i, 'j':j, 'k':k, 'l':l, 'func':improper['func'], 'phase':improper['phase'], 'kd':improper['kd']})


		return newBondList, newPairList, newAngleList, newProperList, newImproperList

			

						
def _bond_updator(resn, bondList, angleList, properList, improperList, amber_flag): # creates newffbonded.itp file

	now = datetime.datetime.now() 
	nowDate = now.strftime('%Y-%m-%d')

	newBonds = []
	newAngles = []
	newPropers = []
	newImpropers = []

	for bond in bondList:
		bonFile = open('ffbonded.itp', 'r')

		flag = 0

		lines = bonFile.readlines()
		bonFile.close()

		infoPresent = False

		for line in lines:
			line = line.strip()
			if len(line) > 0:
				segments = line.split()
				if segments[0] != ';':
					if segments[0] == '[':
						if segments[1] == 'bondtypes':
							flag = 1
							continue
						elif segments[1] == 'angletypes':
							flag = 2
							continue
					if flag == 1:
						if bond['i'] == segments[0] and bond['j'] == segments[1]:
							infoPresent = True							
							break
						if bond['i'] == segments[1] and bond['j'] == segments[0]:
							infoPresent = True							
							break
					if flag > 1:
						break

		if not infoPresent and bond not in newBonds:
			newBonds.append(bond)


	for angle in angleList:
		bonFile = open('ffbonded.itp', 'r')

		flag = 0

		lines = bonFile.readlines()
		bonFile.close()

		infoPresent = False

		for line in lines:
			line = line.strip()
			if len(line) > 0:
				segments = line.split()
				if segments[0] != ';':
					if segments[0] == '[':
						if segments[1] == 'angletypes':
							flag = 1
							continue
						elif segments[1] == 'dihedraltypes':
							flag = 2
							continue
					if flag == 1:
						if angle['i'] == segments[0] and angle['j'] == segments[1] and angle['k'] == segments[2]:
							infoPresent = True							
							break
						if angle['k'] == segments[0] and angle['j'] == segments[1] and angle['i'] == segments[2]:
							infoPresent = True							
							break
					if flag > 1:
						break

		if not infoPresent and angle not in newAngles:
			newAngles.append(angle)


	if amber_flag:
		for proper in properList:
			bonFile = open('ffbonded.itp', 'r')
			flag = 0	
			lines = bonFile.readlines()
			bonFile.close()
	
			infoPresent = False

			for line in lines:
				line = line.strip()
				if len(line) > 0:
					segments = line.split()
					if segments[0] != ';':
						if segments[0] == '[':
							if segments[1] == 'dihedraltypes' and flag == 0:
								flag = 1
								continue
							elif segments[1] == 'dihedraltypes' and flag == 1:
								flag = 2
								continue
						if flag == 2:
							if proper['i']==segments[0] and proper['j']==segments[1] and proper['k']==segments[2] and proper['l']==segments[3]:
								infoPresent = True							
								break
							if proper['i']==segments[3] and proper['j']==segments[2] and proper['k']==segments[1] and proper['l']==segments[0]:
								infoPresent = True							
								break
							
			if not infoPresent and proper not in newPropers:
				newPropers.append(proper)

		for improper in improperList:
			bonFile = open('ffbonded.itp', 'r')
			flag = 0	
			lines = bonFile.readlines()
			bonFile.close()
	
			infoPresent = False

			for line in lines:
				line = line.strip()
				if len(line) > 0:
					segments = line.split()
					if segments[0] != ';':
						if segments[0] == '[':
							if segments[1] == 'dihedraltypes' and flag == 0:
								flag = 1
								continue
							elif segments[1] == 'dihedraltypes' and flag == 1:
								flag = 2
								continue
						if flag == 1 and improper['k'] == segments[2]:
							iList = [improper['i'], improper['j'], improper['l']]
							jList = [segments[0], segments[1], segments[3]]
							found = 0
							for i in iList:
								if i in jList:
									found += 1
							if found == 3:
								infoPresent = True							
								break
							
			if not infoPresent and improper not in newImpropers:
				newImpropers.append(improper)

	else:

		for proper in properList:
			bonFile = open('ffbonded.itp', 'r')
			flag = 0	
			lines = bonFile.readlines()
			bonFile.close()
	
			infoPresent = False

			for line in lines:
				line = line.strip()
				if len(line) > 0:
					segments = line.split()
					if segments[0] != ';':
						if segments[0] == '[':
							if segments[1] == 'dihedraltypes' and flag == 0:
								flag = 1
								continue
							elif segments[1] == 'dihedraltypes' and flag == 1:
								flag = 2
								continue
						if flag == 1:
							if proper['i']==segments[0] and proper['j']==segments[1] and proper['k']==segments[2] and proper['l']==segments[3]:
								infoPresent = True							
								break
							if proper['i']==segments[3] and proper['j']==segments[2] and proper['k']==segments[1] and proper['l']==segments[0]:
								infoPresent = True							
								break
							
			if not infoPresent and proper not in newPropers:
				newPropers.append(proper)

		for improper in improperList:
			bonFile = open('ffbonded.itp', 'r')
			flag = 0	
			lines = bonFile.readlines()
			bonFile.close()
	
			infoPresent = False

			for line in lines:
				line = line.strip()
				if len(line) > 0:
					segments = line.split()
					if segments[0] != ';':
						if segments[0] == '[':
							if segments[1] == 'dihedraltypes' and flag == 0:
								flag = 1
								continue
							elif segments[1] == 'dihedraltypes' and flag == 1:
								flag = 2
								continue
						if flag == 2 and improper['i'] == segments[0]:
							iList = [improper['j'], improper['k'], improper['l']]
							jList = [segments[1], segments[2], segments[3]]
							found = 0
							for i in iList:
								if i in jList:
									found += 1
							if found == 3:
								infoPresent = True							
								break
							
			if not infoPresent and improper not in newImpropers:
				newImpropers.append(improper)


	thingsToWrite = []

	bonFile = open('ffbonded.itp', 'r')

	lines = bonFile.readlines()

	if len(newBonds) > 0:
		newBondCounter = True
	else:
		newBondCounter = False
	if len(newAngles) > 0:
		newAngleCounter = True
	else:
		newAngleCounter = False
	if len(newPropers) > 0:
		newProperCounter = True
	else:
		newProperCounter = False
	if len(newImpropers) > 0:	
		newImproperCounter = True
	else:
		newImproperCounter = False

	flag = 0

	for line in lines:
		if len(line) > 0:
			if not line.startswith(';'):
				if line.startswith('['):
					segments = line.split()
					if segments[1] == 'bondtypes':
						flag = 1
						thingsToWrite.append(line)
						continue
					elif segments[1] == 'angletypes':
						flag = 2
						thingsToWrite.append(line)
						continue
					elif segments[1] == 'dihedraltypes':
						if flag == 2:
							flag = 3
							thingsToWrite.append(line)
							continue
						else:
							flag = 4
							thingsToWrite.append(line)
							continue
		
		if flag == 1 and newBondCounter:
			thingsToWrite.append('\n; bond parameters from '+resn+' ('+nowDate+')\n')

			for newbond in newBonds:
				if amber_flag:
					thingsToWrite.append(f"  {newbond['i']:<3}{newbond['j']:<3}{newbond['func']:>9}{newbond['b0']:>11.5f}{newbond['kb']:>11.1f}\n")
				else:
					thingsToWrite.append(f"{newbond['i']:>8}{newbond['j']:>9}{newbond['func']:>6}{newbond['b0']:>13.8f}{newbond['kb']:>13.2f}\n")
			newBondCounter = False

		elif flag == 2 and newAngleCounter:
			thingsToWrite.append('\n; angle parameters from '+resn+' ('+nowDate+')\n')

			for newangle in newAngles:
				if amber_flag:
					thingsToWrite.append(f"{newangle['i']:<4}{newangle['j']:<4}{newangle['k']:<4}{newangle['func']:>10}{newangle['th0']:>10.3f}{newangle['cth']:>12.3f}\n")
				else:
					thingsToWrite.append(f"{newangle['i']:>8}{newangle['j']:>9}{newangle['k']:>9}{newangle['func']:>6}{newangle['th0']:>13.6f}{newangle['cth']:>12.6f}{newangle['S0']:>13.8f}{newangle['Kub']:>13.2f}\n")

			newAngleCounter = False

		elif amber_flag and flag == 3 and newImproperCounter:
			thingsToWrite.append('\n; improper dihedral parameters from '+resn+' ('+nowDate+')\n')

			for newimproper in newImpropers:
				thingsToWrite.append(f"{newimproper['i']:<4}{newimproper['j']:<4}{newimproper['k']:<4}{newimproper['l']:<4}{newimproper['func']:>6}{newimproper['phase']:>12.2f}{newimproper['kd']:>12.5f}{newimproper['pn']:>6}\n")

			newImproperCounter = False
			
		elif amber_flag and flag == 4 and newProperCounter:
			thingsToWrite.append('\n; proper dihedral parameters from '+resn+' ('+nowDate+')\n')
			
			for newproper in newPropers:
				thingsToWrite.append(f" {newproper['i']:<4}{newproper['j']:<4}{newproper['k']:<4}{newproper['l']:<4}{newproper['func']:>3}{newproper['phase']:>12.1f}{newproper['kd']:>13.5f}{newproper['pn']:>6}\n")
			newProperCounter = False

		elif not amber_flag and flag == 3 and newProperCounter:
			thingsToWrite.append('\n; proper dihedral parameters from '+resn+' ('+nowDate+')\n')
			for newproper in newPropers:
				thingsToWrite.append(f"{newproper['i']:>8}{newproper['j']:>9}{newproper['k']:>9}{newproper['l']:>9}{newproper['func']:>6}{newproper['phase']:>13.6f}{newproper['kd']:>13.6f}{newproper['pn']:>6}\n")
			newProperCounter = False

		elif not amber_flag and flag == 4 and newImproperCounter:
			thingsToWrite.append('\n; improper dihedral parameters from '+resn+' ('+nowDate+')\n')

			for newimproper in newImpropers:
				thingsToWrite.append(f"{newimproper['i']:>8}{newimproper['j']:>9}{newimproper['k']:>9}{newimproper['l']:>9}{newimproper['func']:>6}{newimproper['phase']:>13.6f}{newimproper['kd']:>13.6f}\n")
			newImproperCounter = False

		else:
			thingsToWrite.append(line)

	if (len(newBonds)==0) and (len(newAngles)==0) and (len(newPropers)==0) and (len(newImpropers)==0):
		return 0
	
	newbonFile = open('newffbonded.itp', 'w')
	
	for line in thingsToWrite:
		newbonFile.write(line)

	newbonFile.close()

	print('Wrote '+str(len(newBonds))+' bond types in newffbonded.itp.')
	print('Wrote '+str(len(newAngles))+' angle types in newffbonded.itp.')
	print('Wrote '+str(len(newPropers))+' proper dihedral types in newffbonded.itp.')
	print('Wrote '+str(len(newImpropers))+' improper types in newffbonded.itp.')

	return 0						


def _nonbond_updator(resn, atomList, pairList, amber_flag): # creates newffnonbonded.itp file
	
	now = datetime.datetime.now() 
	nowDate = now.strftime('%Y-%m-%d')

	newAtoms = []
	newPairs = []

	for atom in atomList:
		nonbonFile = open('ffnonbonded.itp', 'r')
	
		flag = 0

		lines = nonbonFile.readlines()
		nonbonFile.close()

		infoPresent = False

		for line in lines:
			line = line.strip()
			if len(line) > 0:
				segments = line.split()
				if segments[0] != ';':
					if segments[0] == '[':
						if segments[1] == 'atomtypes':
							flag = 1
							continue
						elif segments[1] == 'pairtypes':
							flag = 2
							continue
					if flag == 1:
						if atom['name'] == segments[0]:
							infoPresent = True							
							break
					if flag > 1:
						break

		if not infoPresent and atom not in newAtoms:
			newAtoms.append(atom)	

		
	for pair in pairList:
		nonbonFile = open('ffnonbonded.itp', 'r')
	
		flag = 0

		lines = nonbonFile.readlines()
		nonbonFile.close()

		infoPresent = False

		for line in lines:
			line = line.strip()
			if len(line) > 0:
				segments = line.split()
				if segments[0] != ';':
					if segments[0] == '[':
						if segments[1] == 'pairtypes':
							flag = 1
							continue
					if flag == 1:
						if pair['i'] == segments[0] and pair['j'] == segments[1]:
							infoPresent = True							
							break

						elif pair['i'] == segments[1] and pair['j'] == segments[0]:
							infoPresent = True

		if not infoPresent and pair not in newPairs:
			newPairs.append(pair)



	thingsToWrite = []

	nonbonFile = open('ffnonbonded.itp', 'r')

	lines = nonbonFile.readlines()

	if len(newAtoms) > 0:
		newAtomCounter = True
	else:
		newAtomCounter = False
	if len(newPairs) > 0:
		newPairCounter = True
	else:
		newPairCounter = False

	for line in lines:
		if len(line) > 0:
			if not line.startswith(';'):
				if line.startswith('['):
					segments = line.split()
					if segments[1] == 'atomtypes':
						flag = 1
						thingsToWrite.append(line)
						continue
					elif segments[1] == 'pairtypes':
						flag = 2
						thingsToWrite.append(line)
						continue

		
		if flag == 1 and newAtomCounter:
			thingsToWrite.append('\n; Atom parameters from '+resn+' ('+nowDate+')\n')

			for newatom in newAtoms:
				if amber_flag:
					thingsToWrite.append(f"{newatom['name']:<2}{newatom['anum']:>12}{newatom['mass']:>12.3f}   0.0000{newatom['ptype']:>3}{newatom['sigma']:>14}{newatom['epsilon']:>13}\n")
				else:
					thingsToWrite.append(f"{newatom['name']:>6}{newatom['anum']:>6}{newatom['mass']:>13.6f}    0.000{newatom['ptype']:>6}{newatom['sigma']:>16.12f}{newatom['epsilon']:>11.7f}\n")
			newAtomCounter = False

		elif flag == 2 and newPairCounter:
			thingsToWrite.append('\n; Pair parameters from '+resn+' ('+nowDate+')\n')
			
			for newpair in newPairs:
				if not amber_flag:
					thingsToWrite.append(f"{newpair['i']:>6}{newpair['j']:>7}{newpair['func']:>6}{newpair['sig']:>16.12f}{newpair['eps']:>16.12f}\n")

			newPairCounter = False

		else:
			thingsToWrite.append(line)		

	if (len(newAtoms) == 0) and (len(newPairs) == 0):
		return 0
	
	newnonbonFile = open('newffnonbonded.itp', 'w')
	
	for line in thingsToWrite:
		newnonbonFile.write(line)

	newnonbonFile.close()

	print('Wrote '+str(len(newAtoms))+' atom types in newffnonbonded.itp.')
	print('Wrote '+str(len(newPairs))+' pair types in newffnonbonded.itp.')

	return 0	


def _atp_updator(resn, atomList): # appends atomtype data to the atomtypes.atp file in case of AMBER forcefield

	now = datetime.datetime.now() 
	nowDate = now.strftime('%Y-%m-%d')

	atpFile = open('atomtypes.atp', 'r')

	lines = atpFile.readlines()
	atpFile.close()

	newAtomType = []

	for atom in atomList:
		infoPresent = False
		for line in lines:
			line = line.strip()
			if len(line) > 0:
				segments = line.split()
				if segments[0] == atom['name']:
					infoPresent = True
		
		if not infoPresent:
			 newAtomType.append(atom)

	if len(newAtomType) == 0:
		return 0

	atpFile = open('atomtypes.atp', 'a')

	atpFile.write('\n; atomtypes parameters from '+resn+' ('+nowDate+')\n')	

	for atom in newAtomType:
		atpFile.write(f"{atom['name']:<10}{atom['mass']:>16.5f}      ; atomtypes from GAFF2\n")

	atpFile.close()

	print("Wrote "+str(len(newAtomType))+' atom types in atomtypes.atp.')

	return 0	
	

def main():
	"""
	main body
	"""
	resn, mol2File, itpFile, prmFile, forcefield, amber_flag = _read_args()

	oldAtoms, atoms, oldBonds, bonds, pairs, angles, propers, impropers, Impropers, atomList, bondList, pairList, angleList, properList, improperList, H_roots = _cap_remover(mol2File, itpFile, prmFile, amber_flag)
	
	_rtp_updator(resn, atoms, bonds, Impropers, forcefield, amber_flag)

	H_group = _hydrogen_sorter(atoms, oldBonds, H_roots, amber_flag)

	_hdb_updator(resn, H_group, forcefield, amber_flag)

	newBondList, newPairList, newAngleList, newProperList, newImproperList = _list_updator(oldAtoms, atoms, pairs, bonds, angles, propers, impropers, bondList, pairList, angleList, properList, improperList, amber_flag)
	
	_bond_updator(resn, newBondList, newAngleList, newProperList, newImproperList, amber_flag)

	_nonbond_updator(resn, atomList, newPairList, amber_flag)

	if amber_flag:
		_atp_updator(resn, atomList)

	print("*******************************************************************************************************************")
	print("All is well. Now add your residue as a protein type in residuetypes.dat in the gromacs/share/gromacs/top directory.")
	print("*******************************************************************************************************************")

	return 0

if __name__ == "__main__":
	main()
