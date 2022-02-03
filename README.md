# addNewResidue.py

https://github.com/kimjc95/addNewResidue.py

This code adds custom-made amino acids to the GROMACS forcefield directory. 

Currently supporting AMBER and CHARMM forcefields.

This code was written by Joo-Chan Kim at Molecular Synthetic Biology Laboratory in Korea Advanced Institute of Science and Technology (http://msbl.kaist.ac.kr).

If you have any questions or improvements to make, please contact me through kimjoochan@kaist.ac.kr .

When you supply the .mol2 file of your amino acid with N-acetyl cap and C-methylamine cap, this code will
1. remove ACE and NME caps
2. compensate charge difference by spreading it over all sidechain atoms
3. relabel hydrogen atoms according to the connectivity
4. add new parameters to the aminoacids.rtp, aminoacids.hdb, and atomtypes.atp files
5. create newffbonded.itp, newffnonbonded.itp files if there is any new parameter to add to ffbonded.itp and ffnonbonded.itp files

*WARNING: This code only works for non-terminal residues! Also, it simply add/subtracts charge change resulting from removing ACE/NME caps.
        
 So the resulting atomic partial charges may not be accurate!
        
 For the accurate parameterization, use other methods such as quantum mechanical ones.

Guide : 

## 1. Modify / Create your residue in the protein .pdb file using structure viewers such as PyMOL.

Name the heavy atoms in your custom-made residue.

The atom names should be unique in the residue and must be less than 4 letters long!

You do not have to add hydrogen atoms.

Then save your .pdb structure with custom residue (without hydrogens).

## 2. Save your residue's structure as MOL2 format.

For simplicity, remove all residues other than i-1, i, i+1 th ones from the .pdb file you made in step 1.

Change (i-1)th residue into acetyl group.

Change (i+1)th residue into N-methylamine group.

You may now add hydrogens according to the pH. You don't have to name them, since addNewResidue.py will rename them afterwards.

(Make sure to check your hydrogens to have appropriate residue index! ((i-1) for ACE hydrogens, i for your residue, (i+1) for NME hydrogens))

Then save as the .mol2 file.

## 3. Optimize your structure using simple molecular mechanics methods.

You may use external programs such as Chem3D or Avogadro.

or you may use the optimize plugin from PyMOL (https://pymolwiki.org/index.php/Optimize)

Save the energy minimized structure as .mol2 file.

## 4. Create parameters for your residue using appropriate tools!

### For AMBER forcefield, use acpype to generate GAFF2 parameters for your residue.

checkout (https://github.com/alanwilter/acpype)

To run the acpype, you have to unify residue names (change ACE and NME's name into your residue's) in your input .mol2 file.

But do not unify the residue indices! addNewResidue.py differentiates the ACE & NME caps from your residue by the residue index info in your .mol2 file!

### For CHARMM forcefield, use CHARMM-GUI to generate CGenFF parameters for your residue.

checkout (https://charmm-gui.org/)

Sometimes CHARMM-GUI may not properly print out the improper dihedral infos for phenyl rings in your residue, so be aware of that!

Place the resulting folder in your working directory.

## 5. In the same directory, download the latest version of AMBER/CHARMM forcefield files from GROMACS website.

checkout (https://www.gromacs.org/Downloads/User_contributions/Force_fields)

Please note that there is an error in amber14sb_OL15.ff_corrected-Na-cation-params.tar.gz file made by mabraham, 08:29, 30 Aug 2019.

you should manually change line in ffbonded.itp file's improper [ dihedraltypes ] section 

        CT  CV  CC  NA       4      180.00     4.60240     2
        
into

        CT  CC  CV  NA       4      180.00     4.60240     2

to avoid errors you will face when running the simulation of systems containing HID residues.

## 6. Pass your .mol2 file and acpype/charmm-gui folder through the addNewResidue.py code.

for AMBER forcefield:

        python addNewResidue.py -m NEW.mol2 -f amber14.ff -n NAME -a NEW.acpype
        
for CHARMM forcefield:
        
        python addNewResidue.py -m NEW.mol2 -f charmm36.ff -n NAME -c charmm-gui-result

For the detailed info about the input format, use -h or --help flag.

If all is well, the aminoacids.rtp, aminoacids.hdb and atomtypes.atp files will not require further changes.

Check the newly created newffbonded.itp and newffnonbonded.itp files.

In case where there is no new bond / nonbonded interaction parameters to add, the newffbonded.itp and newffnonbonded.itp files will not be created.

If everything seems fine, change their names into ffbonded.itp and ffnonbonded.itp.

## 7. Add your residue's name to the /gromacs/share/gromacs/top/residuetypes.dat, and set the type as Protein.

## 8. Process your .pdb file from step 1 through gmx pdb2gmx.

Enjoy simulation!
 
