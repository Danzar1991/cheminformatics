from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw



reagent = Chem.MolFromSmiles("CC(C)[C@H](O)[C@@H](CC)C(=O)[C@@H](C)C(C)(C)CCCC(C)C(C)C")
check_pattern = Chem.MolFromSmiles("CC(C)C(=O)C(C)CO")
add_pattern_list = ["CC(Cl)C(=O)C(C)CO", "CC(Br)C(=O)C(C)CO", "OCC(C)C(=O)C(C)CO"]
if not reagent.HasSubstructMatch(check_pattern):
    print("Молекула не может быть получена альдольно-кротоновой конденсацией")
elif any([reagent.HasSubstructMatch(Chem.MolFromSmiles(pattern)) for pattern in add_pattern_list]):
    print("Молекула не может быть получена альдольно-кротоновой конденсацией")
else:
    def part_volume(molecule):
        tert_butyl_moiety = Chem.MolFromSmiles("CC(C)(C)F")
        ter = Chem.AddHs(tert_butyl_moiety)
        AllChem.EmbedMolecule(ter)
        volume_tert_butyl = Chem.AllChem.ComputeMolVolume(ter)

        iso_propyl_moiety = Chem.MolFromSmiles("CC(C)F")
        iso = Chem.AddHs(iso_propyl_moiety)
        AllChem.EmbedMolecule(iso)
        volume_iso_propyl = Chem.AllChem.ComputeMolVolume(iso)

        ethyl_moiety = Chem.MolFromSmiles("CCF")
        ethyl = Chem.AddHs(ethyl_moiety)
        AllChem.EmbedMolecule(ethyl)
        volume_ethyl = Chem.AllChem.ComputeMolVolume(ethyl)

        methyl_moiety = Chem.MolFromSmiles("CF")
        methyl = Chem.AddHs(methyl_moiety)
        AllChem.EmbedMolecule(methyl)
        volume_methyl = Chem.AllChem.ComputeMolVolume(methyl)
        if molecule.HasSubstructMatch(tert_butyl_moiety) and any([molecule.GetAtomWithIdx(index).GetAtomMapNum() for index
                                                              in molecule.GetSubstructMatch(tert_butyl_moiety)]):
            return volume_tert_butyl
        elif molecule.HasSubstructMatch(iso_propyl_moiety) and any([molecule.GetAtomWithIdx(index).GetAtomMapNum() for index
                                                              in molecule.GetSubstructMatch(iso_propyl_moiety)]):
            return volume_iso_propyl
        elif molecule.HasSubstructMatch(ethyl_moiety) and any([molecule.GetAtomWithIdx(index).GetAtomMapNum() for index
                                                              in molecule.GetSubstructMatch(ethyl_moiety)]):
            return volume_ethyl
        elif molecule.HasSubstructMatch(methyl_moiety) and any([molecule.GetAtomWithIdx(index).GetAtomMapNum() for index
                                                              in molecule.GetSubstructMatch(methyl_moiety)]):
            return volume_methyl


    def MolWithoutIsotopes(mol):
        atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]
        for atom, isotope in atom_data:
            if isotope:
                atom.SetIsotope(0)
        return mol

    rxn1 = AllChem.ReactionFromSmarts('[OH1][C:7]([*:1])[C:8]([*:2])[C:6][C:5]([*:3])([*:4])'
                                      '>>[OH1][C:7]([111*:1])[C:8]([222*:2])[C:6][C:5]([333*:3])([444*:4])')
    product = rxn1.RunReactants((reagent, ))
    Metka = Chem.MolFromSmiles(AllChem.MolToSmiles(product[0][0]))

    rxn2 = AllChem.ReactionFromSmarts('[OH1][C:7]([111*:1])[C:8]([222*:2])[C:6][C:5]([333*:3])([444*:4])'
                                      '>>[OH1][C:7]([111*:1])[C:8]([222*:2])[C:6][C:5](F)([444*:4]).F[333*:3]')
    product1 = rxn2.RunReactants((Metka, ))
    R3 = Chem.MolFromSmiles(AllChem.MolToSmiles(product1[0][1]))
    for atom in R3.GetAtoms():
        if atom.GetAtomicNum() == 9:
            for a_nei in atom.GetNeighbors():
                a_nei.SetAtomMapNum(9)
    R3_H = Chem.AddHs(R3)
    AllChem.EmbedMolecule(R3_H)
    V3 = part_volume(R3_H)

    rxn3 = AllChem.ReactionFromSmarts('[OH1][C:7]([111*:1])[C:8]([222*:2])[C:6][C:5]([333*:3])([444*:4])'
                                      '>>[OH1][C:7]([111*:1])[C:8]([222*:2])[C:6][C:5]([333*:3])(F).F[444*:4]')
    product1 = rxn3.RunReactants((Metka, ))

    R4 = Chem.MolFromSmiles(AllChem.MolToSmiles(product1[0][1]))
    for atom in R4.GetAtoms():
        if atom.GetAtomicNum() == 9:
            for a_nei in atom.GetNeighbors():
                a_nei.SetAtomMapNum(9)

    R4_H = Chem.AddHs(R4)
    AllChem.EmbedMolecule(R4_H)
    V4 = part_volume(R4_H)




    pattern1 = Chem.MolFromSmarts('C[C@H](O)[C@@H](C)C(=O)[C@@H]C')
    pattern2 = Chem.MolFromSmarts('C[C@@H](O)[C@H](C)C(=O)[C@@H]C')

    if reagent.HasSubstructMatch(pattern1, useChirality=True) and (V3 < V4):
            rxnFinal = AllChem.ReactionFromSmarts(
                '[OH1][C@:7]([111*:1])[C@:8]([222*:2])[C:6][C@:5]([333*:3])([444*:4])>>'
                'O=[12*:7]([12*:1]).[C@:5]([12*:3])([12*:4])[12*:6][12*:8][12*:2]')
            products = rxnFinal.RunReactants((Metka,))
            Draw.ShowMol(MolWithoutIsotopes(products[0][0]))
            Draw.ShowMol(MolWithoutIsotopes(products[0][1]))

    if reagent.HasSubstructMatch(pattern2, useChirality=True) and (V3 > V4):
            rxnFinal = AllChem.ReactionFromSmarts(
                '[OH1][C@:7]([111*:1])[C@:8]([222*:2])[C:6][C@@:5]([333*:3])([444*:4])>>'
                'O=[12*:7]([12*:1]).[C@@:5]([12*:3])([12*:4])[12*:6][12*:8][12*:2]')
            products = rxnFinal.RunReactants((Metka,))
            Draw.ShowMol(MolWithoutIsotopes(products[0][0]))
            Draw.ShowMol(MolWithoutIsotopes(products[0][1]))
