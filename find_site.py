from Bio import SeqIO

human_site_Zn1 = ["C20", "H76"] # ["C10", "H56"]
human_site_Zn2 = ["H108", "H110", "E112"] # ["H88", "H90", "E92"]
human_site_Zn3 = ["H51", "D94", "E92"] # ["H31", "D74", "E72"]

def get_ref_indices(id, seq, site): # pegar os índices dos resíduos de referência de acordo com a posição deles no alinhamento/arquivo
    slash_i = id.find('/') # a barra indica a seperação entre o nome da proteina e o intervalo da sequência
    if (slash_i == -1): return -1 # caso não tenha a barra no id
    
    seq_interval = id[(slash_i + 1):]
    hyphen_i = seq_interval.find('-')
    seq_i = int(seq_interval[:hyphen_i]) # número do primeiro resíduo na sequência
    site_found = []

    # conferir se o número dos resíduos em "site" é igual ao número na sequência "seq"
    for i, residue in enumerate(seq):
        if (len(site_found) == len(site)): break # caso já tenha encontrado todo o sitio na sequência "seq"
        if (residue != '-'):
            if (residue + str(seq_i) in site): 
                res_align_num = i # índice do resíduo no arquivo
                site_found.append(residue + str(res_align_num))

            seq_i += 1
    return site_found    

def find_ref_residues(seq, residues_ref_indices):
    for residue_ref in residues_ref_indices: # procurar cada resíduo de referência na sequência "seq"
        ref_pos = int(residue_ref[1:]) # posição de referência do resíduo
        if (residue_ref[0] != seq[ref_pos]): return -1 # basta um dos resíduos não estar presente para que não tenha o sítio de ligação na sequência
    
    return 0

proteins_found = [] # lista para armazenar as proteínas encontradas que possuem sitios de zinco
for i, record in enumerate(SeqIO.parse("PF00576_panther.fasta", "fasta")):
    record_id = record.id
    record_seq = record.seq

    if (i == 0):
        residues_ref_indices = get_ref_indices(record_id, record_seq, human_site_Zn1)
        continue

    found = find_ref_residues(record_seq, residues_ref_indices)
    if (found == 0): proteins_found.append(record_id)

for protein in proteins_found: print(protein)