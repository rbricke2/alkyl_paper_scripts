"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Makes .rtp files of nucleotides containing an alkyl-phosphate: guanine, cytosine, adenine, and thymine with decyl- or ethyl-phosphate.
    
    Assumes the following files are in the working directory:
        * DP0.mol2 : input CGenFF file of dimethyl decyl-phosphate (modification located on non-bridging oxygen O2P)
        * DP0.rtp  : output CGenFF file of dimethyl decyl-phosphate (modification located on non-bridging oxygen O2P);
                     file is located in the GROMACS format folder
        * DP1.mol2 : input CGenFF file of dimethyl decyl-phosphate (modification located on non-bridging oxygen O1P)
        * DP1.rtp  : output CGenFF file of dimethyl decyl-phosphate (modification located on non-bridging oxygen O1P);
                     file is located in the GROMACS format folder
        * EP0.mol2 : input CGenFF file of dimethyl ethyl-phosphate (modification located on non-bridging oxygen O2P)
        * EP0.rtp  : output CGenFF file of dimethyl ethyl-phosphate (modification located on non-bridging oxygen O2P);
                     file is located in the GROMACS format folder
        * EP1.mol2 : input CGenFF file of dimethyl ethyl-phosphate (modification located on non-bridging oxygen O1P)
        * EP1.rtp  : output CGenFF file of dimethyl ethyl-phosphate (modification located on non-bridging oxygen O1P);
                     file is located in the GROMACS format folder
        * DTN.mol2 : input CGenFF file of deoxythymidine
        * DTN.rtp  : output CGenFF file of deoxythymidine; file is located in the GROMACS format folder
        * DCN.mol2 : input CGenFF file of deoxycytidine
        * DCN.rtp  : output CGenFF file of deoxycytidine; file is located in the GROMACS format folder
        * DAN.mol2 : input CGenFF file of deoxyadenosine
        * DAN.rtp  : output CGenFF file of deoxyadenosine; file is located in the GROMACS format folder
        * DGN.mol2 : input CGenFF file of deoxyguanosine
        * DGN.rtp  : output CGenFF file of deoxyguanosine; file is located in the GROMACS format folder
    
    Outputted files:
        * AD0.rtp  : adenine with decyl-phosphate (modification located on non-bridging oxygen O2P)
        * AD1.rtp  : adenine with decyl-phosphate (modification located on non-bridging oxygen O1P)
        * TD0.rtp  : thymine with decyl-phosphate (modification located on non-bridging oxygen O2P)
        * TD1.rtp  : thymine with decyl-phosphate (modification located on non-bridging oxygen O1P)
        * CD0.rtp  : cytosine with decyl-phosphate (modification located on non-bridging oxygen O2P)
        * CD1.rtp  : cytosine with decyl-phosphate (modification located on non-bridging oxygen O1P)
        * GD0.rtp  : guanine with decyl-phosphate (modification located on non-bridging oxygen O2P)
        * GD1.rtp  : guanine with decyl-phosphate (modification located on non-bridging oxygen O1P)
        * AE0.rtp  : adenine with ethyl-phosphate (modification located on non-bridging oxygen O2P)
        * AE1.rtp  : adenine with ethyl-phosphate (modification located on non-bridging oxygen O1P)
        * TE0.rtp  : thymine with ethyl-phosphate (modification located on non-bridging oxygen O2P)
        * TE1.rtp  : thymine with ethyl-phosphate (modification located on non-bridging oxygen O1P)
        * CE0.rtp  : cytosine with ethyl-phosphate (modification located on non-bridging oxygen O2P)
        * CE1.rtp  : cytosine with ethyl-phosphate (modification located on non-bridging oxygen O1P)
        * GE0.rtp  : guanine with ethyl-phosphate (modification located on non-bridging oxygen O2P)
        * GE1.rtp  : guanine with ethyl-phosphate (modification located on non-bridging oxygen O1P)
"""

################################################################################################
#
# FUNCTIONS
#
################################################################################################

def get_atom_names(mol2_file):
    atom_names       = []
    store_atom_names = False
    with open(mol2_file, "r") as f:
        for line in f:
            if 'BOND' in line:
                break
            elif store_atom_names:
                atom_names.append(line.split()[1])
            elif 'ATOM' in line:
                store_atom_names = True

    return atom_names

def list_to_rtp_atom_line(line_split):
    return str(line_split[0].rjust(9)+line_split[1].rjust(7)+line_split[2].rjust(9)+line_split[3].rjust(3)+'\n')

def change_atom_names_and_neutralize(rtp_file, atom_names):
    new_rtp         = []
    atoms           = False
    bonds           = False
    impropers       = False
    i               = 0
    atom_names_dict = {}
    with open(rtp_file, "r") as f:
        for line in f:
            new_line   = line
            line_split = line.split()
            if impropers:
                if line_split: # make sure list isn't empty
                    for j in range(4): # impropers require four atoms
                        line_split[j] = atom_names_dict[line_split[j]] # update atom name
                    new_line = str(line_split[0].rjust(9)+line_split[1].rjust(7)+line_split[2].rjust(7)+line_split[3].rjust(7)+'\n')
            elif '[ impropers ]' in line:
                bonds     = False
                impropers = True
            elif bonds:
                for j in range(2): # bonds require two atoms
                    line_split[j] = atom_names_dict[line_split[j]] # update atom name
                new_line = str(line_split[0].rjust(9)+line_split[1].rjust(7)+'\n')
            elif '[ bonds ]' in line:
                atoms = False
                bonds = True
            elif atoms:
                atom_names_dict[line_split[0]] = atom_names[i]
                line_split[0]                  = atom_names[i]

                # adjust charge of the atoms O5', C5', O3', and C3' so that the total charge
                # of the nucleotide is zero
                if atom_names[i] == "O5'":
                    line_split[2] = str(float(line_split[2]) + 0.004)
                elif atom_names[i] == "O3'":
                    line_split[2] = str(float(line_split[2]) + 0.003)
                elif atom_names[i] == "C5'":
                    line_split[2] = str(float(line_split[2]) + 0.003)
                elif atom_names[i] == "C3'":
                    line_split[2] = str(float(line_split[2]) + 0.003)

                new_line = list_to_rtp_atom_line(line_split)
                i       += 1
            elif '[ atoms ]' in line:
                atoms = True

            new_rtp.append(str(new_line))

    return new_rtp

def change_atom_types(new_rtp_list, base):
    # Define dictionaries for atom types
    common = {
        "C5'" : "CN8B",
        "H5'1": "HN8",
        "H5'2": "HN8",
        "C4'" : "CN7",
        "H4'" : "HN7",
        "O4'" : "ON6",
        "C1'" : "CN7B",
        "H1'" : "HN7",
        "C2'" : "CN8",
        "H2'1": "HN8",
        "H2'2": "HN8",
        "C3'" : "CN7",
        "H3'" : "HN7",
        "O3'" : "ON2"
    }
    
    thymine = {
        "N1" : "NN2B",
        "C2" : "CN1T",
        "O2" : "ON1",
        "N3" : "NN2U",
        "H3" : "HN2",
        "C4" : "CN1",
        "O4" : "ON1",
        "C5" : "CN3T",
        "C5M": "CN9",
        "C7" : "CN9",
        "H51": "HN9",
        "H52": "HN9",
        "H53": "HN9",
        "H71": "HN9",
        "H72": "HN9",
        "H73": "HN9",
        "C6" : "CN3",
        "H6" : "HN3"
    }
    
    adenine = {
        "N9" : "NN2",
        "C5" : "CN5",
        "N7" : "NN4",
        "C8" : "CN4",
        "H8" : "HN3",
        "N1" : "NN3A",
        "C2" : "CN4",
        "H2" : "HN3",
        "N3" : "NN3A",
        "C4" : "CN5",
        "C6" : "CN2",
        "N6" : "NN1",
        "H61": "HN1",
        "H62": "HN1"
    }
    
    cytosine = {
        "N1" : "NN2",
        "C2" : "CN1",
        "O2" : "ON1C",
        "N3" : "NN3",
        "C4" : "CN2",
        "N4" : "NN1",
        "H41": "HN1",
        "H42": "HN1",
        "C5" : "CN3",
        "H5" : "HN3",
        "C6" : "CN3",
        "H6" : "HN3"
    }
    
    guanine = {
        "N9" : "NN2B",
        "C4" : "CN5",
        "N2" : "NN1",
        "H21": "HN1",
        "H22": "HN1",
        "N3" : "NN3G",
        "C2" : "CN2",
        "N1" : "NN2G",
        "H1" : "HN2",
        "C6" : "CN1",
        "O6" : "ON1",
        "C5" : "CN5G",
        "N7" : "NN4",
        "C8" : "CN4",
        "H8" : "HN3"
    }
    
    if base == "T": # thymine
        merged = common | thymine
    elif base == "C": # cytosine
        merged = common | cytosine
    elif base == "G": # guanine
        merged = common | guanine
    else: # adenine
        merged = common | adenine
    
    new_rtp           = []
    change_atom_types = False
    for line in new_rtp_list:
        new_line = line
        if '[ bonds ]' in line:
            change_atom_types = False
        if change_atom_types:
            line_split = line.split()
            atom_name  = line_split[0]
            if atom_name in merged:
                line_split[1] = merged[atom_name]
            new_line = list_to_rtp_atom_line(line_split)
        elif '[ atoms ]' in line:
            change_atom_types = True
        
        new_rtp.append(str(new_line))
    
    return new_rtp

def comment_out_capping_atoms(list_rtp_file):
    new_rtp = []
    for line in list_rtp_file:
        new_line   = line
        line_split = line.split()
        
        if len(line_split) >= 2:
            if 'R' in line_split[0] or 'R' in line_split[1]:
                # removing capping atoms
                new_line = ''.join((';', new_line))
            elif line_split[0] == "O3'" and line_split[1] == "P":
                # removing the bond between the phosphorus atom and the
                # 3’ methoxy oxygen of the dimethyl alkyl-phosphate
                new_line = ''.join((';', new_line))
        
        new_rtp.append(str(new_line))

    return new_rtp

def get_atoms_and_bonds(list_rtp_file):
    atom_line = False
    bond_line = False
    atoms     = []
    bonds     = []
    for line in list_rtp_file:
        if '[ impropers ]' in line:
            break
        if bond_line:
            bonds.append(str(line))
        elif '[ bonds ]' in line:
            bond_line = True
            atom_line = False
        if atom_line:
            atoms.append(str(line))
        elif '[ atoms ]' in line:
            atom_line = True
        
    return atoms, bonds

def get_residue_name(list_rtp_file):
    residue_name = ""
    for i in range(len(list_rtp_file)):
        if '[ atoms ]' in list_rtp_file[i]:
            residue_name = str(list_rtp_file[i-1].split()[1])
            break
        
    return residue_name

def merge_rtp_files(new_deoxyribonucleoside_rtp, new_dimethyl_alkyl_phosphate_rtp):
    merged_rtp = []
    
    backbone_residue = get_residue_name(new_dimethyl_alkyl_phosphate_rtp)
    base_residue     = get_residue_name(new_deoxyribonucleoside_rtp)
    
    if 'E' in backbone_residue:
        dimethyl_alkyl_phosphate_desc = "dimethyl ethyl-phosphate"
        modification                  = "an ethyl-phosphate group"
    else:
        dimethyl_alkyl_phosphate_desc = "dimethyl decyl-phosphate"
        modification                  = "a decyl-phosphate group"
    
    if base_residue == "DTN":
        deoxyribonucleoside_desc = "deoxythymidine"
        base                     = "thymine"
        new_residue              = "T"+backbone_residue[0]+backbone_residue[-1]
    elif base_residue == "DCN":
        deoxyribonucleoside_desc = "deoxycytidine"
        base                     = "cytosine"
        new_residue              = "C"+backbone_residue[0]+backbone_residue[-1]
    elif base_residue == "DGN":
        deoxyribonucleoside_desc = "deoxyguanosine"
        base                     = "guanine"
        new_residue              = "G"+backbone_residue[0]+backbone_residue[-1]
    else:
        deoxyribonucleoside_desc = "deoxyadenosine"
        base                     = "adenine"
        new_residue              = "A"+backbone_residue[0]+backbone_residue[-1]


    additional_header_lines = [str("; Merged .rtp files of " + dimethyl_alkyl_phosphate_desc + " and " + deoxyribonucleoside_desc + " to\n"),
                               str("; create central " + base + " fragment containing " + modification + ".\n"),
                               str("; Created using make_rtp_file.py, written by Rachel Bricker.\n\n")]
    new_deoxyribonucleoside_rtp = additional_header_lines + new_deoxyribonucleoside_rtp
    
    backbone_atoms, backbone_bonds = get_atoms_and_bonds(new_dimethyl_alkyl_phosphate_rtp)
    
    additional_bonds = [str("; ATOM CONNECTIVITIES BETWEEN THE METHOXY OXYGENS OF THE\n"),
                        str("; " + dimethyl_alkyl_phosphate_desc.upper() + " AND THE C5' AND C3' ATOMS OF\n"),
                        str("; THE " + deoxyribonucleoside_desc.upper() + "\n"),
                        str("C3'".rjust(9)+"O3'".rjust(7)+'\n'),
                        str("C5'".rjust(9)+"O5'".rjust(7)+'\n'),
                        str("; BACKBONE CONNECTION TO THE NEXT RESIDUE\n"),
                        str("O3'".rjust(9)+"+P".rjust(7)+'\n')]
    
    for i in range(len(new_deoxyribonucleoside_rtp)):
        line       = new_deoxyribonucleoside_rtp[i]
        new_line   = line
        line_split = line.split()
        if line_split:
            if ('[ ' + base_residue + ' ]') in line:
                new_line = str("[ " + new_residue + " ]\n; CAPPING ATOMS CONTAIN THE LETTER 'R'\n")
            elif '[ atoms ]' in line:
                new_line = new_line + str("; ATOMS FROM " + deoxyribonucleoside_desc.upper() + "\n")
            elif '[ bonds ]' in line:
                charge_group = 0
                for j in range(len(backbone_atoms)):
                    backbone_line_split = backbone_atoms[j].split()

                    if backbone_line_split[-1] == "0":
                        charge_group = int(new_deoxyribonucleoside_rtp[i-1].split()[-1])+1

                    backbone_line_split[-1] = str(charge_group)
                    
                    if len(backbone_line_split) == 5:
                        backbone_atoms[j] = backbone_line_split[0] + list_to_rtp_atom_line(backbone_line_split[1:])
                    else:
                        backbone_atoms[j] = list_to_rtp_atom_line(backbone_line_split)
                    charge_group += 1

                backbone_atoms.insert(0, str("; ATOMS FROM " + dimethyl_alkyl_phosphate_desc.upper() + " GROUP\n"))
                backbone_bonds.insert(0, str("; BONDS FROM " + dimethyl_alkyl_phosphate_desc.upper() + " GROUP\n"))
                    
                merged_rtp = merged_rtp + backbone_atoms
                new_line   = new_line + str("; BONDS FROM " + deoxyribonucleoside_desc.upper() + "\n")
            elif '[ impropers ]' in line:
                merged_rtp = merged_rtp + backbone_bonds + additional_bonds
                new_line   = new_line + str("; IMPROPERS FROM " + deoxyribonucleoside_desc.upper() + "\n")
        merged_rtp.append(str(new_line))
    
    return merged_rtp, new_residue

def prepare_rtp_file_for_merging(rtp_file, mol2_file, base='T'):
    new_rtp = change_atom_names_and_neutralize(rtp_file, get_atom_names(mol2_file))
    new_rtp = change_atom_types(new_rtp, base)
    new_rtp = comment_out_capping_atoms(new_rtp)
    
    return new_rtp

def get_total_charge(merged_rtp_list):
    atom_line = False
    charge    = 0
    for line in merged_rtp_list:
        if '[ bonds ]' in line:
            break
        if atom_line:
            line_split = line.split()
            if len(line_split) == 4:
                charge    += float(line_split[2])
        elif '[ atoms ]' in line:
            atom_line = True
    
    return charge

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    dimethyl_alkyl_phosphates = ["DP0", "DP1", "EP0", "EP1"]
    deoxyribonucleosides      = {"DTN" : "T",
                                 "DCN" : "C",
                                 "DAN" : "A",
                                 "DGN" : "G"}
    for i in range(len(dimethyl_alkyl_phosphates)):
        for residue in deoxyribonucleosides:
            # CGenFF input (.mol2 files)
            dimethyl_alkyl_phosphate_mol2    = dimethyl_alkyl_phosphates[i] + ".mol2"
            deoxyribonucleoside_mol2         = residue + ".mol2"

            # CGenFF output (.rtp files)
            dimethyl_alkyl_phosphate_rtp     = dimethyl_alkyl_phosphates[i] + ".rtp"
            deoxyribonucleoside_rtp          = residue + ".rtp"
            
            # modify .rtp files so that they are ready for merging
            new_deoxyribonucleoside_rtp      = prepare_rtp_file_for_merging(deoxyribonucleoside_rtp, deoxyribonucleoside_mol2, deoxyribonucleosides[residue])
            new_dimethyl_alkyl_phosphate_rtp = prepare_rtp_file_for_merging(dimethyl_alkyl_phosphate_rtp, dimethyl_alkyl_phosphate_mol2)

            # merge .rtp files
            merged_rtp, new_residue = merge_rtp_files(new_deoxyribonucleoside_rtp, new_dimethyl_alkyl_phosphate_rtp)
            
            # create new .rtp file for residue
            with open(new_residue + ".rtp", "w+") as rtp:
                for line in merged_rtp:
                    rtp.write(line)
                print("File " + new_residue + ".rtp was created!\n")
                # print the sum of charges (should be zero)
                print("The total charge of the nucleotide is: " + str(get_total_charge(merged_rtp)) + "\n")

if __name__ == "__main__": 
    main()