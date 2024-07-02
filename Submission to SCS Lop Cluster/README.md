# Submission to SCS Lop Cluster
<h2 align="center">
  
  This documentation is to detail how to submit Gaussian Jobs on the SCS Lop Cluster
  <br>
  
  📄📝➡️💻🧬📈👨‍💻
</h2>

<div>
  
**'.xyz'** files are simple text-based files that represent the atomic coordinates for a molecule. Programs that output '.xyz' files are those involved in exploring molecular modeling such as CREST (Conformer-Rotamer Ensemble Sampling Tool).

**'.gjf'** files are used to run quantum chemistry calculations on Gaussian software.  <br> <br>

Sample Code is provided at **XYZ_to_GJF.ipynb** in the *xyz File to gjf File* folder or click here [XYZ_to_GJF](https://github.com/SelvinTo/CompChem-Resources/blob/574c3a8b56e183f148e992afca3a438864fb602b/xyz%20File%20to%20gjf%20File/XYZ_to_GJF.ipynb)

<br>
<br>

<details>
  <summary> Click to view Step by Step Instructions </summary>
  
  ## Instructions
  
  1. **Step 1**: Copy the code as provided in XYZ_to_GJF.ipynb, this will convert a single '.xyz' file to a '.gjf' file.  

    import os

    def crest_xyz_to_gjf(xyz_file, output_folder, method='B3LYP', basis_set='6-31G(d)', charge=0, multiplicity=1):   #Can change Method/Basis_Set here
      
      """
      Convert a CREST XYZ file with multiple conformers to separate Gaussian input (GJF) files,
      and store them in a specified output folder.

      Parameters:
      xyz_file (str or file-like object): Path to the input CREST XYZ file or a file-like object.
      output_folder (str): Path to the folder where the output GJF files will be saved.
      method (str): The computational method to use in Gaussian (default: B3LYP).
      basis_set (str): The basis set to use in Gaussian (default: 6-31G(d)).
      charge (int): The total charge of the molecule (default: 0).
      multiplicity (int): The multiplicity of the molecule (default: 1).
      """
     
      # Create the output folder if it doesn't exist
      if not os.path.exists(output_folder):
          os.makedirs(output_folder)

      # Open the XYZ file
      if isinstance(xyz_file, str):
          xyz_file = open(xyz_file, 'r')

      lines = xyz_file.readlines()
      xyz_file.close()

      index = 0
      num_lines = len(lines)
      conformer_index = 0

      while index < num_lines:
          # Read number of atoms
          num_atoms = int(lines[index].strip())
          index += 1

          # Read energy or comments
          energy_or_comment = lines[index].strip()
          index += 1

          # Extract atomic coordinates for the conformer
          atoms = []
          for _ in range(num_atoms):
              parts = lines[index].split()
              element = parts[0]
              x, y, z = map(float, parts[1:4])
              atoms.append((element, x, y, z))
              index += 1

          # Write the GJF file for the current conformer
          gjf_filename = os.path.join(output_folder, f"conformer_{conformer_index + 1}.gjf")
          with open(gjf_filename, 'w') as gjf_file:
              gjf_file.write('%mem=16GB\n')
              gjf_file.write('%nprocshared=16\n')
              gjf_file.write(f'# opt freq {method}/{basis_set} empiricaldispersion=gd3bj integral=ultrafine \n\n')    #Gaussian Input_Line/Parameters 
              gjf_file.write(f'{energy_or_comment}\n\n')
              gjf_file.write(f'{charge} {multiplicity}\n')
              for atom in atoms:
                  element, x, y, z = atom
                  gjf_file.write(f'{element:>2}{x:>28.8f}{y:>14.8f}{z:>14.8f}\n')
              gjf_file.write('\n')

          conformer_index += 1

  2. **Step 2**: Adjust any Gaussian Input parameters as neccesary. These two lines control the Method, Basis Set, and gaussian input.
     
       ***Geometry Optimization*** can be perfomed at the PBE-D3(BJ)/6-31+G(d,p) level of theory. PBE-D3(BJ) being the method and 6-31+G(d,p) being the basis_set.

       ***Single Point Energy Calculations*** can be perfomed at the PBE0-D3(BJ)/def2-TZVP level of theory. PBE0-D3(BJ) being the method and def2-TZVP being the basis_set.

       Additional terms that can be appended at the end depending on usage included (***nmr=giao, prop=efg, pop=nbo, freq=noraman***)
     
     
    def crest_xyz_to_gjf(xyz_file, output_folder, method='B3LYP', basis_set='6-31G(d)', charge=0, multiplicity=1):   #Can change Method/Basis_Set here

    gjf_file.write(f'# opt freq {method}/{basis_set} empiricaldispersion=gd3bj integral=ultrafine \n\n')    #Gaussian Input_Line/Parameters 
    
  
  4. **Step 3**: Example Usage 

    xyz_file_path = 'Test.xyz'   #Put path to xyz file here
    output_folder = 'Test.gjf'   #Put path to output folder here 
    crest_xyz_to_gjf(xyz_file_path, output_folder)
  
</details>
