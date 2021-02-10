
Instructions to generate frames of array assembly which can be used to create cool videos as shown in:
https://twitter.com/Arielbs100/status/1346931500336365570


1. Create conda environment:
  a. Add repository to conda config file: 
      conda config --add channels schrodinger
  b. Create a new environment:
      conda env create -f pymolVid.yml
  c. Activate the environemt:
      conda activate pymolVid37
 
2. Updating pymol for production:
  a. Open the file bashCMD.txt:
      vi bashCMD.txt
  b. Update line 2 with the correct path to the relevant pymol and project
  c. Add alias to .bashrc:
      head -3 bashCMD.txt >> ~/.bashrc
  d. source ~/.bashrc

3. run example (this can be time consuming process in which it will dump .png files by default to outputSlides/test/ unless specified otherwise ):
   pmlVid Acomp.pdb Bcomp.pdb 

4. Create .mp4 file from a folder containing the enumerated file created in step 3 (path to input slides and output video file are specified in source/slides2videoRev.py [lines 22, 23]):
    python source/slides2videoRev.py   
   To see the result open  output/Vid/test.mp4  
      

Note:
  This project is purposed to aid beginers in simulation based mulecular visualization, please contact with any bugs in running the example. 

Related work:
   This video was purposed to visualize mulecular assembly of de novo designed binary protein arrays (https://www.nature.com/articles/s41586-020-03120-8) 

   
