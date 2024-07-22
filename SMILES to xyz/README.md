# XYZ File to GJF File
<h2 align="center">
  
  This documentation is to detail how to transform SMILES into '.xyz' files 
  <br>
  
  📝 ➡️ 📐 ➡️ 📄
</h2>

<div>
  
**'.xyz'** files are simple text-based files that represent the atomic coordinates for a molecule. Programs that output '.xyz' files are those involved in exploring molecular modeling such as CREST (Conformer-Rotamer Ensemble Sampling Tool).

**SMILES** strings are "Simplified Molecular Input Line Entry System" that are used to translate a chemical's three-dimensional structure into a string of symbols that is easily understood by computer software.  <br> <br>

Sample Code is provided at **SMILES_to_xyz.ipynb** in the *SMILES to xyz* folder or click here [SMILES_to_xyz](https://github.com/SelvinTo/CompChem-Resources/blob/15850a64462448708b81373abd147e1497e90566/SMILES%20to%20xyz/SMILES_to_xyz.ipynb)

<br>
<br>

<details>
  <summary> Click to view Step by Step Instructions </summary>
  
  ## Instructions
  
  1. **Step 1**: Copy the code as provided in SMILES_to_xyz.ipynb, this will convert a single SMILES string into a single '.xyz' file 

<div> 
   
    pip install rdkit
  
</div>
   
    from rdkit import Chem 

<div> 
   
    from rdkit.Chem import AllChem
  
</div>
   
    def smiles_to_xyz(smiles, output_file):
    mol = Chem.MolFromSmiles(smiles)
    mol_h = Chem.AddHs(mol)  # Adding Hydrogens
    AllChem.EmbedMolecule(mol_h, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    AllChem.MMFFOptimizeMolecule(mol_h) # Computing 3D coordinates
    
    with open(output_file, 'w') as f:
        f.write(f'{mol_h.GetNumAtoms()}\n')
        f.write(f'Generated from SMILES: {smiles}\n')
        conf = mol_h.GetConformer()
        for i in range(mol_h.GetNumAtoms()):
            atom = mol_h.GetAtomWithIdx(i)
            symbol = atom.GetSymbol()
            x, y, z = conf.GetAtomPosition(i)
            f.write(f'{symbol:>2} {x:>18.8f} {y:14.8f} {z:>14.8f}\n')     
    
  2. **Step 2**: Input SMILE string to convert 
<div> 
   
    smiles_RuPhos = 'CC(C)OC(C=CC=C1OC(C)C)=C1C(C=CC=C2)=C2P(C3CCCCC3)C4CCCCC4'  # Example SMILES string (Aspirin)
    output_file = 'RuPhos.xyz'

    smiles_to_xyz(smiles_RuPhos, output_file)
  
</div>
The output xyz file should look like this:  
    <img src="RuPhos.xyz" width="90%"/>  <br>

  
  4. **Step 3**: Example Usage 

    xyz_file_path = 'Test.xyz'   #Put path to xyz file here
    output_folder = 'Test.gjf'   #Put path to output folder here 
    crest_xyz_to_gjf(xyz_file_path, output_folder)
  
</details>
