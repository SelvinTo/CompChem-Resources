# Submission to SCS Lop Cluster
<h2 align="center">
  
  This documentation is to detail how to submit Gaussian Jobs on the SCS Lop Cluster
  <br>
  
  📄📝➡️💻🧬📈👨‍💻
</h2>

<div>
  
The **SCS Lop Cluster** is available for use by anyone in the School of Chemical Sciences at UIUC. To request access and learn more information, you can vist the official site here: https://answers.uillinois.edu/scs/page.php?id=104365

<br>
This documentation is to help detail how to run Gaussian Jobs for both Windows and Mac Users.
<br>
<br>

<details>
  <summary> Click to view Windows Users Step by Step Instructions </summary>
  
  ## Instructions
  
  1. **Step 1**: Uploading a folder from local computer to your cluster folder.
     
     a. Open Windows Powershell (type powershell in search bar) <br>
     b. In the command line, type the following and press enter: scp -r [filelocation] netID@lop.scs.illinois.edu:/home/NetID <br>
     C. You should be prompted with your netID password. Enter that into the command line. You should then see all the files successfully copied, as shown below <br>
       
    

  3. **Step 2**: Adjust any Gaussian Input parameters as neccesary. These two lines control the Method, Basis Set, and gaussian input.
     
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
<div>
<details>
<summary> Click to view Mac Users Step by Step Instructions </summary>

  ## Instructions
  
 4. **Step 3**: Example Usage 

    xyz_file_path = 'Test.xyz'   #Put path to xyz file here
    output_folder = 'Test.gjf'   #Put path to output folder here 
    crest_xyz_to_gjf(xyz_file_path, output_folder)
