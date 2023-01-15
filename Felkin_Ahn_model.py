from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw



def MolWithoutIsotopes(mol):
    atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]
    for atom, isotope in atom_data:
        if isotope:
            atom.SetIsotope(0)
    return mol

# Целевая молекула
reagent = Chem.MolFromSmiles('O[C@@](C)(C(C)(C)C)[C@@](C(C)C)(C1=CC=CC=C1)C2C(C)C(C)C3=CC=CC=C3C2')
Draw.ShowMol(reagent)

rxn1 = AllChem.ReactionFromSmarts('[OH1][C:6]([*:1])([*:2])[C:7]([*:3])([*:4])([*:5])>>[OH1][C:6]([111*:1])([222*:2])[C:7]([333*:3])([444*:4])([555*:5])')
product = rxn1.RunReactants((reagent, ))
Metka = Chem.MolFromSmiles(AllChem.MolToSmiles(product[0][0]))
#Metka - Целевая молекула с "изотопными метками"

# Вспомогательная реакция №2.
# Механически отрывает заместитель R3 от целевой молекулы с метками -- с меткой "333".
# Места разрыва заменяет атомами фтора. Добавляет к молекуле R3-F водороды, генерирует конформацию R3-F, рассчитывает объем R3-F -- V3.

rxn2 = AllChem.ReactionFromSmarts('[OH1][C:6]([111*:1])([222*:2])[C:7]([333*:3])([444*:4])([555*:5])>>[OH1][C:6]([111*:1])([222*:2])[C:7](F)([444*:4])([555*:5]).F[333*:3]')
product1 = rxn2.RunReactants((Metka,))
R3 = Chem.MolFromSmiles(AllChem.MolToSmiles(product1[0][1]))
R3_H = Chem.AddHs(R3)
AllChem.EmbedMolecule(R3_H)
V3 = Chem.AllChem.ComputeMolVolume(R3_H)

# Вспомогательная реакция №3 -- то же с заместителем R4.

rxn3 = AllChem.ReactionFromSmarts('[OH1][C:6]([111*:1])([222*:2])[C:7]([333*:3])([444*:4])([555*:5])>>[OH1][C:6]([111*:1])([222*:2])[C:7]([333*:3])(F)([555*:5]).F[444*:4]')
product1 = rxn3.RunReactants((Metka,))
R4 = Chem.MolFromSmiles(AllChem.MolToSmiles(product1[0][1]))
R4_H = Chem.AddHs(R4)
AllChem.EmbedMolecule(R4_H)
V4 = Chem.AllChem.ComputeMolVolume(R4_H)

# Вспомогательная реакция №4 -- то же с заместителем R5.

rxn4 = AllChem.ReactionFromSmarts('[OH1][C:6]([111*:1])([222*:2])[C:7]([333*:3])([444*:4])([555*:5])>>[OH1][C:6]([111*:1])([222*:2])[C:7]([333*:3])([444*:4])(F).F[555*:5]')
product1 = rxn4.RunReactants((Metka,))
R5 = Chem.MolFromSmiles(AllChem.MolToSmiles(product1[0][1]))
R5_H = Chem.AddHs(R5)
AllChem.EmbedMolecule(R5_H)
V5 = Chem.AllChem.ComputeMolVolume(R5_H)

# Все возможные стереопаттерные, которые могут заматчится на целевую молекулу с метками (см презентацию).

pattern1 = Chem.MolFromSmarts('[OH1][C@:6]([111*:1])([222*:2])[C@:7]([333*:3])([444*:4])([555*:5])')
pattern2 = Chem.MolFromSmarts('[OH1][C@@:6]([111*:1])([222*:2])[C@:7]([333*:3])([444*:4])([555*:5])')
pattern3 = Chem.MolFromSmarts('[OH1][C@:6]([111*:1])([222*:2])[C@@:7]([333*:3])([444*:4])([555*:5])')
pattern4 = Chem.MolFromSmarts('[OH1][C@@:6]([111*:1])([222*:2])[C@@:7]([333*:3])([444*:4])([555*:5])')

# Применение правила Фелкина-Ана.
# Осуществление финальной стереселективной реакции в зависимости от отношения объемов заместителей.

if Metka.HasSubstructMatch(pattern1, useChirality=True) or Metka.HasSubstructMatch(pattern4, useChirality=True):
    if (V3 < V4 < V5) or (V4 < V5 < V3) or (V5 < V3 < V4):
        rxnFinal = AllChem.ReactionFromSmarts('[OH1][C@:6]([111*:1])([222*:2])[C@:7]([333*:3])([444*:4])([555*:5])>>O=[C:6]([12*:1])[C@:7]([12*:3])([12*:4])[12*:5].I[#12][12*:2]')
        products = rxnFinal.RunReactants((Metka,))
        Draw.ShowMol(MolWithoutIsotopes(products[0][0]))
        Draw.ShowMol(MolWithoutIsotopes(products[0][1]))
    if (V3 < V5 < V4) or (V4 < V3 < V5) or (V5 < V4 < V3):
        rxnFinal = AllChem.ReactionFromSmarts('[OH1][C@:6]([111*:1])([222*:2])[C@:7]([333*:3])([444*:4])([555*:5])>>O=[C:6]([12*:2])[C@:7]([12*:3])([12*:4])[12*:5].I[#12][12*:1]')
        products = rxnFinal.RunReactants((Metka,))
        Draw.ShowMol(MolWithoutIsotopes(products[0][0]))
        Draw.ShowMol(MolWithoutIsotopes(products[0][1]))

if Metka.HasSubstructMatch(pattern2, useChirality=True) or Metka.HasSubstructMatch(pattern3, useChirality=True):
    if  (V3 < V4 < V5) or (V4 < V5 < V3) or (V5 < V3 < V4):
        rxnFinal = AllChem.ReactionFromSmarts('[OH1][C@:6]([111*:1])([222*:2])[C@:7]([333*:3])([444*:4])([555*:5])>>O=[C:6]([12*:2])[C@:7]([12*:3])([12*:4])[12*:5].I[#12][12*:1]')
        products = rxnFinal.RunReactants((Metka,))
        Draw.ShowMol(MolWithoutIsotopes(products[0][0]))
        Draw.ShowMol(MolWithoutIsotopes(products[0][1]))
    if (V3 < V5 < V4) or (V4 < V3 < V5) or (V5 < V4 < V3):
        rxnFinal = AllChem.ReactionFromSmarts('[OH1][C@:6]([111*:1])([222*:2])[C@:7]([333*:3])([444*:4])([555*:5])>>O=[C:6]([12*:1])[C@:7]([12*:3])([12*:4])[12*:5].I[#12][12*:2]')
        products = rxnFinal.RunReactants((Metka,))
        Draw.ShowMol(MolWithoutIsotopes(products[0][0]))
        Draw.ShowMol(MolWithoutIsotopes(products[0][1]))


