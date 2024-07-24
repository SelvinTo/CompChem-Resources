# UMAP
<h2 align="center">
  
  This documentation details how to perform a UMAP Projection
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
  

  3. **Step 3**: Input SMILES 

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
  <img src="Screenshot 2024-07-24 115240.png" width="60%"/>  <br> 

4. **Step 4**: Create UMAP Projection  

</div>
   
     # Select features for UMAP projection (ECFP6 fingerprints)
    features = df_ecfp6.drop(columns=['smiles'])
    
<div> 
   
    from umap import UMAP
  
</div>

    # Perform UMAP projection
    umap = UMAP(n_components=2)  # Set the number of dimensions for projection
    umap_projection = umap.fit_transform(features)

</div>
   
     # Visualize the data
     plt.scatter(umap_projection[:, 0], umap_projection[:, 1], alpha=0.5)
     plt.title('UMAP Projection of ECFP6 Fingerprints')
     plt.xlabel('UMAP Component 1')
     plt.ylabel('UMAP Component 2')
     plt.show()
  <img src="Screenshot 2024-07-24 115335.png" width="60%"/>  <br> 
<div> 
   
    # Create a color map based on the indices you want to highlight
    indices_to_highlight = [1,100, 802, 903, 104, 105, 306, 207, 308, 409, 501, 111, 212, 613, 814, 415, 616, 117, 118]  # Specify the indices you want to highlight
    colors = ['red' if i in indices_to_highlight else 'blue' for i in range(len(df))]
  
</div>   
  
    # Perform UMAP projection
    umap_model = UMAP(n_components=2)
    umap_projection = umap_model.fit_transform(features)

<div> 

    # Visualize the data with the color map
    plt.figure(figsize=(8, 6))
    plt.scatter(umap_projection[:, 0], umap_projection[:, 1], c=colors, alpha=0.5)
    plt.title('UMAP Projection of ECFP6 Fingerprints')
    plt.xlabel('UMAP Component 1')
    plt.ylabel('UMAP Component 2')
    plt.show()
 <img src="Screenshot 2024-07-24 115350.png" width="60%"/>  <br> 

</h2>

</div>  
</details>
