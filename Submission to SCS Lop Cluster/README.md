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
     
     **a.** Open Windows Powershell (type powershell in search bar) <br><br>
     **b.** In the command line, type the following and press enter: *scp -r [local-filelocation] netID@lop.scs.illinois.edu:/home/NetID* <br><br>
     **c.** You should be prompted with your netID password. Enter that into the command line. You should then see all the files successfully copied, as shown below <br>
     <img src="Screenshot 2024-07-22 113459.png" width="90%"/>    
     Hint:<br>
- To copy just a single file and not a directory/folder, remove '-r' from the command line.<br>
- To copy the local-filelocation of your folder/file, you can right click on the file and click "copy as path" <br>

2. **Step 2**: Connect to the cluster
     
      **a.** In the command line, type " *ssh -Y netID@lop.scs.illinois.edu* " <br><br>
      **b.** You will then be prompted to enter your Illinois password <br><br>
      **c.** If you see this window, that means you have successfully connected to the cluster!<br>
 
      <img src="Screenshot 2024-07-22 113608.png" width="90%"/>  <br>
  
3. **Step 3**: Submitting Jobs to the Cluster
     
      **a.** In the command line, type *ls* and press enter.This will show you all the files in your directory. You should see the folders that you uploaded to the cluster in the previous step  <br>
      <img src="Screenshot 2024-07-22 113632.png" width="90%"/>  

     **b.** Type *cd [foldername]* to enter the folder you just uploaded. Type *ls* and press enter to see all the files within that folder. (This is not neccessary if only a single file was uploaded instead of a directory/folder.   <br>
      <img src="Screenshot 2024-07-22 113658.png" width="90%"/>  
 
     **c.** In the command line, type *module load gaussian/g16* and press enter<br><br>
     Congrats! You are now ready to run gaussian job files (.gjf) <br><br>
     **d.** In the command line, type *submit-g16 -n 16 -q [clustername] [filepath]* and press enter<br>
     - This command tells the computer to submit a guassian16 calculation with 16 cores to this processor, and the job file is found at this location <br>
     - More information about the clustername can be found in the Cluster Queue Section; the filepath is the location of the file on your cluster directory. <br>
       <img src="Screenshot 2024-07-22 113722.png" width="90%"/> <br>
       
     To see if the job was submitted successfully, type *ls*. You should now see .log files for each job ran. You can repeat the command line for each .gjf file you wish to submit a job for. Hint: you can press the up arrow on the keyboard to load the code you entered previously, which you can edit with the arrow keys/backspace. This makes it easier than typing the same line every time 😊    <br>
      **e.** Congrats! You have submitted gaussian job files on the cluster! <br>
   - To view the status of each job, type *qstat* and press enter <br>
   - To view how much computing time you are using, type *qquota* and press enter <br>
4. **Step 4**: Exiting the Cluster <br>
      **a.** Type *exit* and press enter to disconnect from the cluster before closing the terminal 

</details>
<div>
<details>
<summary> Click to view Mac Users Step by Step Instructions </summary>

  ## Instructions
  
 1. **Step 1**: Uploading a folder from local computer to your cluster folder.
     
     **a.** Open a new windonw in the terminal <br><br>
     **b.** In the command line, type the following and press enter: *scp -r [local-filelocation] netID@lop.scs.illinois.edu:/home/NetID* <br><br>
     **c.** You should be prompted with your netID password. Enter that into the command line. You should then see all the files successfully copied, as shown below <br>
     <img src="Screenshot 2024-07-22 113459.png" width="90%"/>    
     Hint:<br>
- To copy just a single file and not a directory/folder, remove '-r' from the command line.<br>
- To copy the file location on your computer, you can “right click” + “options” the folder and click “Copy File as Pathname”. You can then paste the pathname into the terminal   <br>

2. **Step 2**: Connect to the cluster
     
      **a.** In the command line, type " *ssh -Y netID@lop.scs.illinois.edu* " <br><br>
      **b.** You will then be prompted to enter your Illinois password <br><br>
      **c.** If you see this window, that means you have successfully connected to the cluster!<br>
 
      <img src="Screenshot 2024-07-22 113608.png" width="90%"/>  <br>
  
3. **Step 3**: Submitting Jobs to the Cluster
     
      **a.** In the command line, type *ls* and press enter.This will show you all the files in your directory. You should see the folders that you uploaded to the cluster in the previous step  <br>
      <img src="Screenshot 2024-07-22 113632.png" width="90%"/>  

     **b.** Type *cd [foldername]* to enter the folder you just uploaded. Type *ls* and press enter to see all the files within that folder. (This is not neccessary if only a single file was uploaded instead of a directory/folder.   <br>
      <img src="Screenshot 2024-07-22 113658.png" width="90%"/>  
 
     **c.** In the command line, type *module load gaussian/g16* and press enter<br><br>
     Congrats! You are now ready to run gaussian job files (.gjf) <br><br>
     **d.** In the command line, type *submit-g16 -n 16 -q [clustername] [filepath]* and press enter<br>
     - This command tells the computer to submit a guassian16 calculation with 16 cores to this processor, and the job file is found at this location <br>
     - More information about the clustername can be found in the Cluster Queue Section; the filepath is the location of the file on your cluster directory. <br>
       <img src="Screenshot 2024-07-22 113722.png" width="90%"/> <br>
       
     To see if the job was submitted successfully, type *ls*. You should now see .log files for each job ran. You can repeat the command line for each .gjf file you wish to submit a job for. Hint: you can press the up arrow on the keyboard to load the code you entered previously, which you can edit with the arrow keys/backspace. This makes it easier than typing the same line every time 😊    <br>
      **e.** Congrats! You have submitted gaussian job files on the cluster! <br>
   - To view the status of each job, type *qstat* and press enter <br>
   - To view how much computing time you are using, type *qquota* and press enter <br>
4. **Step 4**: Exiting the Cluster <br>
      **a.** Type *exit* and press enter to disconnect from the cluster before closing the terminal <br>

</details>
<div>
<details>
  <summary> Click to view Cluster Queue Information </summary>
</div>
