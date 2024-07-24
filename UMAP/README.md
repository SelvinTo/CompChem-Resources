# UMAP
<h2 align="center">
  
  This documentation is to detail how to do a UMAP Projection
  <br>
  
  📝 ➡️ 📐 ➡️ 📄
</h2>

<div>
  
**UMAP** projections are a helpful way of visualizing data and is a type of dimensionality reduction. <br> <br>

Sample Code is provided at **UMAP.ipynb** in the *UMAP* folder or click here [UMAP](https://github.com/SelvinTo/CompChem-Resources/blob/ac8fb2cefa1ec46270b0eea8b3e1e198def3ef9f/UMAP/UMAP.ipynb)
<br>
<br>

<details>
  <summary> Click to view Step by Step Instructions </summary>
  
  ## Instructions
  
  1. **Step 1**: Copy the code as provided in UMAP.ipynb, this will create a UMAP projection of a list of SMILES. <br> <br> Install neccessary libraries.  

<div> 
   
    pip install pandas
  
</div>
   
    pip install rdkit

<div> 
   
    pip install umap-learn
  
</div>
    
    pip install matplotlib

<div> 
    
    pip install openpyxl
  
</div>

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap,LinearSegmentedColormap
    from sklearn.cluster import KMeans

    from rdkit import Chem,Geometry
    from rdkit.Chem import rdmolfiles, AllChem, rdMolAlign,rdmolops, Descriptors, Draw

    from collections import Counter
    from rdkit.Chem import AllChem
    from rdkit import Chem, DataStructs
    import numpy as np
    import pandas as pd    
    
  2. **Step 2**: Create ECFP Calculator 
<div> 
   
    #This calculates the molecular vectors from smiles strings
    class ECFPCalculator:
      def __init__(self, smiles):
          self.mols = [Chem.MolFromSmiles(i) for i in smiles]
          self.smiles = smiles

      def mol2fp(self, mol, radius = 3):
          fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius = radius)
          array = np.zeros((1,))
          Chem.DataStructs.ConvertToNumpyArray(fp, array)
          return array

      def compute_ECFP6(self, name):
          bit_headers = ['bit' + str(i) for i in range(2048)]
          arr = np.empty((0,2048), int).astype(int)
          for i in self.mols:
              fp = self.mol2fp(i)
              arr = np.vstack((arr, fp))
          df_ecfp6 = pd.DataFrame(np.asarray(arr).astype(int),columns=bit_headers)
          df_ecfp6.insert(loc=0, column='smiles', value=self.smiles)
          df_ecfp6.to_csv(name[:-4]+'_ECFP6.csv', index=False)
  

  4. **Step 3**: Input SMILES 

<div> 
   
    # Read the Excel file containing SMILES strings
    df = pd.read_excel('UMAP_Test_Smiles.xlsx')

    #Looks for a column named 'smiles' and converts it into a list
    smiles_list = df['smiles'].tolist()

    # Create an instance of the class with the SMILES strings
    ecfp_calculator = ECFPCalculator(smiles_list)

    # Specify the name of the output CSV file
    output_csv_name = "ECFP_test_output.csv"

    # Compute ECFP6 fingerprints and save the results to a CSV file
    ecfp_calculator.compute_ECFP6(output_csv_name)

    # Read the CSV file containing the ECFP6 fingerprints
    df_ecfp6 = pd.read_csv("ECFP_test_output_ECFP6.csv")

    # Verify the contents of the DataFrame
    print(df_ecfp6.head())  # Display the first few rows of the DataFrame
  
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
