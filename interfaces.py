from Bio.PDB import PDBParser, Polypeptide
import numpy

def similarity(interface1, interface2):
    if len(interface1) >= len(interface2):
        greater = len(interface1)
    else:
        greater = len(interface2)

    equal_contacts = 0
    for contact1 in interface1: # contact = [residue1, atom1, residue2, atom2, distance]
        res_name1a = contact1[0].get_resname()
        num1a = contact1[0].id[1]
        res_name1b = contact1[2].get_resname()
        num1b = contact1[2].id[1]

        for contact2 in interface2:
            res_name2a = contact2[0].get_resname()
            num2a = contact2[0].id[1]
            res_name2b = contact2[2].get_resname()
            num2b = contact2[2].id[1]

            if ((res_name1a, num1a) == (res_name2a, num2a) and (res_name1b, num1b) == (res_name2b, num2b)) or ((res_name1a, num1a) == (res_name2b, num2b) and (res_name1b, num1b) == (res_name2a, num2a)):
                equal_contacts += 1
                break
    
    return equal_contacts, greater

def print_similarity(interface1, interface2, equal_contacts, greater):
    model1a = interface1[0][0].get_parent().get_parent()
    full_id1a = interface1[0][1].get_full_id()
    model1b = interface1[0][2].get_parent().get_parent()
    full_id1b = interface1[0][3].get_full_id()

    model2a = interface2[0][0].get_parent().get_parent()
    full_id2a = interface2[0][1].get_full_id()
    model2b = interface2[0][2].get_parent().get_parent()
    full_id2b = interface2[0][3].get_full_id()

    print(f"Similaridade entre as interfaces {full_id1a[2]}{model1a.id}/{full_id1b[2]}{model1b.id} e {full_id2a[2]}{model2a.id}/{full_id2b[2]}{model2b.id}: {equal_contacts/greater} ({equal_contacts}/{greater})")

def get_similarity(interfaces):
    qtd_interfaces = len(interfaces)

    for i, interface in enumerate(interfaces): # itera sobre todos os indices (i) e elementos (interface) de interfaces
        if i == (qtd_interfaces-1):  # última interface
            break

        n = i + 1  # índice da próxima interface        
        while n < qtd_interfaces:
            next_interface = interfaces[n]
            equal_contacts, greater = similarity(interface, next_interface)
            print_similarity(interface, next_interface, equal_contacts, greater)
            n += 1 # índice da próxima interface (na sequência em interfaces)


def compare_chain_distances(chain1, chain2): # compara as distâncias entre os átomos de duas cadeias
    pairs = []  # array para guardar os pares de residuos já encontrados, para evitar repetições na saída
    contacts = []

    for atom1 in chain1:
        residue1 = atom1.get_parent()
        for atom2 in chain2:
            residue2 = atom2.get_parent()
            distance = numpy.linalg.norm(atom2.coord - atom1.coord)

            if distance < 3.5:
                if (residue1, residue2) not in pairs and (residue2, residue1) not in pairs: # verifica se o par de resíduos ainda não foi registrado
                    pairs.append((residue1, residue2))
                    contacts.append([residue1, atom1, residue2, atom2, distance])
    
    return contacts

def print_contacts(interfaces): # imprime os contatos de todas as interfaces
    print('{:^8} {:^8} {:^8} {:^8} | {:^8} {:^8} {:^8} {:^8} | {:^8}'.format('Modelo', 'Cadeia', 'Resíduo', 'Número', 'Modelo', 'Cadeia', 'Resíduo', 'Número', 'Distância'))
    print('-'*90)
    
    for contacts in interfaces:
        for contact in contacts: # contact = [residue1, atom1, residue2, atom2, distance]
            model1 = contact[0].get_parent().get_parent()
            res_name1 = contact[0].get_resname()
            full_id1 = contact[1].get_full_id()

            model2 = contact[2].get_parent().get_parent()
            res_name2 = contact[2].get_resname()
            full_id2 = contact[3].get_full_id()

            distance = contact[4]

            print('{:^8} {:^8} {:^8} {:^8} | {:^8} {:^8} {:^8} {:^8} | {:^8.3f}'.format(model1.id, full_id1[2], res_name1, contact[0].id[1],
                                                                                        model2.id, full_id2[2], res_name2, contact[2].id[1], distance))
        print('Número de Contatos: ', len(contacts))
        print('-'*90)

def compare_chains(chains):
    interfaces = []  # array para armazenar todos os contatos encontrados
    qtd_chains = len(chains)

    for i, chain in enumerate(chains): # itera sobre todos os indices (i) e elementos (chain) de chains
        if i == (qtd_chains-1):  # última cadeia
            break

        n = i + 1  # índice da próxima cadeia
        while n < qtd_chains:
            next_chain = chains[n]
            interface = compare_chain_distances(chain, next_chain)
            interfaces.append(interface)
            n += 1 # índice da próxima cadeia (na sequência em chains)

    #print_contacts(interfaces)
    return interfaces


# carregar a estrutura e extrair as cadeias
parser = PDBParser()
structure = parser.get_structure('Protein', '../PDBs e Fastas/UA/8uw4.pdb')

chains = []
for model in structure:
    for chain in model:
        atoms = []
        for residue in chain:
            for atom in residue:
                if Polypeptide.is_aa(residue.get_resname()):  # considera apenas os resíduos de aminoácidos
                    atoms.append(atom)
        chains.append(atoms)

contacts = compare_chains(chains)
get_similarity(contacts)