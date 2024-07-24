# UMAP
<h2 align="center">
  
  This documentation is to detail how to do a UMAP Projection
  <br>
  
  📝 ➡️ 📐 ➡️ 📄
</h2>

<div>
  
**'.xyz'** files are simple text-based files that represent the atomic coordinates for a molecule. Programs that output '.xyz' files are those involved in exploring molecular modeling such as CREST (Conformer-Rotamer Ensemble Sampling Tool).

**SMILES** strings are "Simplified Molecular Input Line Entry System" that are used to translate a chemical's three-dimensional structure into a string of symbols that is easily understood by computer software.  <br> <br>

Sample Code is provided at **UMAP.ipynb** in the *UMAP* folder or click here [SMILES_to_xyz](https://github.com/SelvinTo/CompChem-Resources/blob/15850a64462448708b81373abd147e1497e90566/SMILES%20to%20xyz/SMILES_to_xyz.ipynb)

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
The output xyz file should look something like this:  
    <img src="Screenshot 2024-07-22 140052.png" width="90%"/>  <br>

  <br>
  
  4. **Step 3**: To convert multiple SMILES into '.xyz' files 

<div> 
   
    pip install pandas
  
</div>
   
    pip install openpyxl
    
<div> 
   
    import pandas as pd
  
</div>
   
    # Read the Excel file
    df = pd.read_excel('example_smiles.xlsx')
    df.head(4)
<img src="Screenshot 2024-07-22 142052.png" width="40%"/>  <br> 
<div> 
   
    for _, row in df.iterrows():
    smile_num = row['Smile #']
    smiles = row['smiles']
    
    # Generate the output file name
    output_file = f'{smile_num}.xyz'
    
    # Convert SMILES to .xyz and save the file
    smiles_to_xyz(smiles, output_file)
<br>
This will have converted all the SMILES string into their respective xyz files.   
<h2 align="center">
   The xyz files can then be opened in software like Avogadro to visualize the molecule or be converted into other file formats. <br><br>
<img src="Screenshot 2024-06-21 at 12.26.24 PM.png" width="60%"/>  <br> 

</h2>

</div>  
</details>
