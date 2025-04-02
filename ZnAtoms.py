from Bio.PDB import PDBParser, Polypeptide
import numpy

def printResidues(residues):
    ids = []
    print('{:^8} {:^8} {:^8} {:<8} {:<10}'.format('Modelo', 'Cadeia', 'Resíduo', 'Número', 'Distância'))
    for residue in residues:
        if residue[3] not in ids:
            ids.append(residue[3])
            print('{:^8} {:^8} {:^8} {:^8} {:^8.3f}'.format(residue[0], residue[1], residue[2], residue[3], residue[4]))

def compareDistances(ZNs, residuesAtms):
    for ZN in ZNs:
        closerResidues = []
        print(f'Átomos próximos ao zinco ', ZN.serial_number, '(', ZN.get_full_id()[2],')')
        #ZN.get_parent().get_parent().get_parent()

        for atm in residuesAtms:
            distanceZN =  numpy.linalg.norm(ZN.coord - atm.coord)

            if distanceZN < 3.5:
                full_id = atm.get_full_id()
                chain = full_id[2]

                residue = atm.get_parent()
                res_name = residue.get_resname()
                res_id = residue.id[1]

                model = residue.get_parent().get_parent()

                resInfo = []
                resInfo.append(model.id)
                resInfo.append(chain)
                resInfo.append(res_name)
                resInfo.append(res_id)
                resInfo.append(distanceZN)

                closerResidues.append(resInfo)
        printResidues(closerResidues)

        print('-'*50, '\n')
            
parser = PDBParser()    
structure = parser.get_structure('Protein', '../PDBs e Fastas/UB/3ssg.pdb1')

residueAtoms = []
ZnAtoms = []

for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.name == 'ZN': 
                   ZnAtoms.append(atom)
                   continue
                elif Polypeptide.is_aa(residue.get_resname()): # considera apenas átomos de resíduos
                    residueAtoms.append(atom)

compareDistances(ZnAtoms, residueAtoms)