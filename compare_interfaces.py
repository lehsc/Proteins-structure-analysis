from Bio.PDB import PDBParser, Polypeptide
from Bio.PDB.Atom import Atom
import numpy

def compareChains(chains):
    print('{:^8} {:^8} {:^8} {:^8} | {:^8} {:^8} {:^8} {:^8} | {:^8}'.format('Modelo','Cadeia', 'Resíduo', 'Número', 'Modelo', 'Cadeia', 'Resíduo', 'Número', 'Distância'))
    print('-'*90)
    
    qtd_chains = len(chains)
    for i, chain in enumerate(chains): # itera sobre todos os indices (i) e elementos (chain) de chains
        if i == (qtd_chains-1): # última cadeia
            break
        
        n = i + 1 # guarda o indice da próxima cadeia
        pairs = [] # array para guardar os pares de residuos já encontrados, para evitar repetições na saída
        while n < qtd_chains:
            next_chain = chains[n]
            qtd_contacts = 0
            for atom1 in chain: # comparar cada átomo da cadeia atual ...
                residue1 = atom1.get_parent()
                for atom2 in next_chain: # ... com cada átomo das cadeias seguintes
                    residue2 = atom2.get_parent()
                    distance = numpy.linalg.norm(atom2.coord - atom1.coord)

                    if distance < 3.5 and [residue1,residue2] not in pairs:
                        qtd_contacts += 1
                        pairs.append([residue1,residue2])

                        full_id1 = atom1.get_full_id()
                        res_name1 = residue1.get_resname()
                        model1 = residue1.get_parent().get_parent()

                        full_id2 = atom2.get_full_id()
                        res_name2 = residue2.get_resname()
                        model2 = residue2.get_parent().get_parent()

                        print('{:^8} {:^8} {:^8} {:^8} | {:^8} {:^8} {:^8} {:^8} | {:^8.3f}'.format(model1.id, full_id1[2], res_name1, residue1.id[1], model2.id, full_id2[2], res_name2, residue2.id[1], distance))
            print('Número de Contatos: ', qtd_contacts)
            print('-'*90)
            n+=1 # indice da próxima cadeia da sequência


parser = PDBParser()
structure = parser.get_structure('Protein', '../PDBs e Fastas/UB/2rox.pdb1')

chains = []
total_residues = 0
for model in structure:
    for chain in model:
        atoms = []
        for residue in chain:
            for atom in residue:
                if Polypeptide.is_aa(residue.get_resname()): #
                    atoms.append(atom)
                    total_residues += 1
        chains.append(atoms)
compareChains(chains)