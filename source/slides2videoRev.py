
 
import sys, os
import numpy as np
from scipy.spatial.transform import Rotation as R
import cv2
import progressbar as pb
from PIL import Image, ImageDraw, ImageFont
from time import strftime
from time import gmtime
import textwrap
cwd = os.getcwd()
p2utilsFolder = os.path.join(cwd,"source")
sys.path.append(p2utilsFolder)
from utils import *



FPS = 30
Reverse = 1 # determine the direction of the enumarated slide files ar sorted  

p2vid = os.path.join(os.getcwd(),"outputVideo/test.mp4")
p2slides = os.path.join(os.getcwd(),"outputSlides/vidSlides")


print("p2vid: ",os.path.basename(p2vid),"\n",p2vid)
imList = sorted([os.path.join(p2slides,ele) for ele in os.listdir(p2slides) if ele[-4:]==".png"], reverse=Reverse)
frameObj = [] 
print(imList)
print("Number of slides: ",len(imList))

[frameObj.append(cv2.imread(p)) for p in imList ]

print("Shape of data: ",type(frameObj),type(frameObj[0]),frameObj[0].shape	 )
pixelsY = frameObj[0].shape[0]
pixelsX = frameObj[0].shape[1]


####totalTime = len(frameObj)/FPS # in seconds 
time=[0,0]
imgY = 150
imgX = 1200
imageSize = [imgX,imgY] ; fontSize = 18 ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[20,100]
text = "Order emerges of chaos - forte of protein engineering"
text = "Protein Engineering Transforming Chaos to Order"
fontSource = '/Library/Fonts/Arial Bold.ttf'
#text2video(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)
imgX = 900
imageSize = [imgX,imgY] ; fontSize = 22 ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[930,90]; time=[0,17]
text = "Ben-Sasson, A. J. et al. \n\nDesign of Biologically Active Binary Protein 2D Materials"
fontSource = '/Library/Fonts/times-new-roman.ttf'
text2video(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)

#time=[17,20]
#text = "\nDesign of Biologically Active Binary Protein 2D Materials"
#fontSource = '/Library/Fonts/times-new-roman.ttf'
#text2video(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)


####
imgY = 150
imgX = 300
imageSize = [imgX,imgY] ; fontSize = 20 ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[80,150] ; time=[17,20]
text = "IPD, UW:\n  Ariel J. Ben-Sasson\n  William Sheffler\n  Justin Decarreau\n  David Baker"
#fontSource = '/Library/Fonts/Arial Bold.ttf'
text2video(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)
####
imageSize = [imgX,imgY] ; fontSize = 20 ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[80,1570]
text = "MRC LMB, Cambridge:\n  Joseph L. Watson\n  Alice Bittleston\n  Emmanuel Derivery"
#fontSource = '/Library/Fonts/Arial Bold.ttf'
text2video(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)
####
imageSize = [imgX,imgY] ; fontSize = 20 ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[280,150]
text = "Biochem., UW:\n  Matthew C. Johnson\n  Justin M. Kollman"
#fontSource = '/Library/Fonts/Arial Bold.ttf'
text2video(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)
####
imageSize = [imgX,imgY] ; fontSize = 20 ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[280,1570]
text = "ISCRM, UW:\n  Logeshwaran\n    Somasundaram\n  Hannele\n    Ruohola-Baker"
#fontSource = '/Library/Fonts/Arial Bold.ttf'
text2video(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)
####
imageSize = [imgX,imgY] ; fontSize = 20 ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[670,150]
text = "PSD, PNNL:\n  Fang Jiao\n  Jiajun Chen\n  James J De Yoreo"
#fontSource = '/Library/Fonts/Arial Bold.ttf'
text2video(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)
####
imageSize = [imgX,imgY] ; fontSize = 20 ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[820,150]
text = "Chem. Eng. Cambridge:\n  Ioanna Mela\n  Clemens F. Kaminski"
#fontSource = '/Library/Fonts/Arial Bold.ttf'
text2video(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)

imageSize = [imgX,imgY] ; fontSize = 20 ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[670,1570]
text = "SIBYLS, LBNL:\n  Greg L. Hura"
#fontSource = '/Library/Fonts/Arial Bold.ttf'
text2video(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)

imageSize = [imgX,imgY] ; fontSize = 20 ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[820,1570]
text = "BCMP, HMS:\n  Andrew A. Drabek\n  Sanchez M. Jarrett\n  Stephen C. Blacklow"
#fontSource = '/Library/Fonts/Arial Bold.ttf'
text2video(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)




def text2video2(frames,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS=30,t=[0,0],Anchor='mm'):
	'''
	This include the opion to change font size with time (an appearing text) 
	Because text become bigger it is important to be able to center the text. 

	'''
	### dimensions 
	pixelsY = frames[0].shape[0]
	pixelsX = frames[0].shape[1]
	yStart = int(position[0]-np.ceil(imageSize[0]/2))
	yEnd = yStart+imageSize[0]
	xStart = int(position[1]-np.ceil(imageSize[1]/2))
	xEnd = xStart+imageSize[1]


	### timing 
	for k in range(len(frames)):
		if (t[0]==0 and t[1]==0):
			frames[k][yStart:yEnd,xStart:xEnd,:] =  np.clip(frames[k][yStart:yEnd,xStart:xEnd,:].astype('int')+imgNP.astype('int'),0,255).astype('uint8')
		elif (k>FPS*t[0] and k<FPS*t[1]):
			ratio = (k-FPS*t[0])/(FPS*(t[1]-t[0]))
			fontSizeR = int((fontSize[1]-fontSize[0])*ratio+fontSize[0])
			img = Image.new('RGB', (imageSize[1], imageSize[0] ), color = bg_color)
			fnt = ImageFont.truetype(fontSource, fontSizeR)
			d = ImageDraw.Draw(img)
			w, h = d.textsize(text, font=fnt)
			d.text((imageSize[0]/2-int(w/2), imageSize[1]/2-int(h/2)), text, fill=textColor, font=fnt,anchor="mb")  # multiline_text
			#d.text((3,10), text, font=fnt, fill=textColor)
			imgNP = np.array(img)

			frames[k][yStart:yEnd,xStart:xEnd,:] =  np.clip(frames[k][yStart:yEnd,xStart:xEnd,:].astype('int')+imgNP.astype('int'),0,255).astype('uint8')
	return frames



### 
imgX=600;imgY=600
imageSize = [imgX,imgY] ; fontSize = [1,40] ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[540,960] ; time=[14,17.5] ; anchor='mm'
text = "Design of Biologically Active \nBinary Protein 2D Materials"
fontSource = '/Library/Fonts/times-new-roman.ttf'
text2video2(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time,anchor)

imageSize = [imgX,imgY] ; fontSize = [40,40] ; bg_color = (0,0,0) ; textColor= (255, 255, 255); position=[540,960] ; time=[17.5,20] ; anchor='mm'
text = "Design of Biologically Active \nBinary Protein 2D Materials"
fontSource = '/Library/Fonts/times-new-roman.ttf'
text2video2(frameObj,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time,anchor)



####imgY = 120
####
####
####
####
####
####
####
####imgY = 30 ; imgX = 130 ; position=[30,20] ; fontSource = '/Library/Fonts/Arial Bold Italic.ttf'
####imageSize = [imgX,imgY]
####clock2video(frameObj,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS=FPS)
####
####
####
####textList=["Design strategy using dihedral (Dn) homooligomers building blocks"\
####		 ,"A 1st component (D3) is positioned at the origin. White arrow indicate components rotation symmetry axes"\
####		 ,"A 2nd component (D2) is positioned such that one of its C2 symmetry axes coincides with one of those of the 1st component"\
####		 ,"D3 (as do D4 and D6) has exactly two unique orientations to meet design constrains"\
####		 ,"D2 is unique and has 6 unique orientations to meet design constrains"\
####		 ,"Docking - components slide into contact. In white: interface regions"\
####		 ,"Interfaces with aligned helices are selected for sequence design"\
####		 ,"The C2 symmetry axis shared by the 2 components and passing through the interface, renders the interface curvature defects resistant"\
####		 ,"Components symmetry guaranties symmetrically repeated interfaces - 3 on each D3 and 2 on each D2 components"\
####		 ,"Propagation to a predefined pattern follows the plane symmetry (p6m) rules"\
####		 ,"We acknowledge the PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC, and python community for the following packages: numpy, scipy, opencv, and pillow"	]
####
####x = 30/FPS
####timeList=[[0,4],[4,8],[8,12],[12,15],[15,19],[19,25],[26,29],[30,35],[36,40],[41,65],[70,82] ]
####timeList = [[l2*x for l2 in  l] for l in timeList]
####position=[70,20] ; fontSize = 16 ; bg_color = (0,0,0) ; textColor= (255, 255, 255)
####
####for text,time in zip(textList,timeList):
####	totalTime = len(frameObj)/FPS # in seconds 
####	#### this part break the line into lines of up to 27 characters#########
####	newList=[]
####	count2 = 0
####	for word,count in  zip(text.split(" "),range(len(text.split(" ")))):
####	    if count==0:
####	        newList.append(word)
####	    elif len(newList[count2]+word)<27:
####	        newList[count2]=newList[count2]+" "+word
####	    else:
####	        newList.append(word)
####	        count2+=1
####	newList2 = "\n".join(newList)
####	########################################################################
####	imgY = 30*len(newList) ; imgX = 250 
####	imageSize = [imgX,imgY] 
####	fontSource = '/Library/Fonts/Arial Bold.ttf'
####	text2video(frameObj,newList2,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS,time)








video = cv2.VideoWriter(p2vid, cv2.VideoWriter_fourcc(*'MP4V'), fps=FPS,frameSize=(pixelsX,pixelsY) ) #frameObj[0].size)
for image in frameObj:
	video.write(image)
video.release()
cv2.destroyAllWindows()



# ./bin/pip install torch-scatter==latest+cpu -f https://pytorch-geometric.com/whl/torch-1.5.0.html
# ./bin/pip install torch-sparse==latest+cpu -f https://pytorch-geometric.com/whl/torch-1.5.0.html
# ./bin/pip install torch-cluster==latest+cpu -f https://pytorch-geometric.com/whl/torch-1.5.0.html
# ./bin/pip install torch-spline-conv==latest+cpu -f https://pytorch-geometric.com/whl/torch-1.5.0.html
# ./bin/pip install torch-geometric




















