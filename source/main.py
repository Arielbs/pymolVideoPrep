

import pymol 
from pymol import cmd
import sys, os
import numpy as np
cwd = os.getcwd()
print("cwd: ",cwd)
p2utilsFolder = os.path.join(cwd,"source")
sys.path.append(p2utilsFolder)
from utils import *


cmd.set('ray_opaque_background', 1)
#### set intial paths objects 

p2source = os.path.join(cwd,"input/pdbs")
p2f = os.path.join(cwd,"outputSlides/test") # path to slide folder
pdpList = sys.argv[3:] # list of input pdb files to be used 
print("Input PDBs list: ",pdpList)
[ cmd.load(os.path.join(p2source,Name)) for Name in pdpList if Name[-4:]==".pdb"]


#### set view properties 
pixelsX = 1920
pixelsY = 1080
view0 = [-0.0,    1.0,    0.0,\
         -1.0,   -0.0,   -0.0,\
         -0.0,   -0.0,    1.0,\
          0.0,    0.0, -1350.0,\
         -0.0,    0.0,    0.0,\
         -14000.0, 14000.0,  -20.0]

viewFinal = [
    -0.0,    1.0,    0.0,\
    -1.0,   -0.0,   -0.0,\
    -0.0,   -0.0,    1.0,\
     0.0,    0.0, -4500.0,\
    -0.0,    0.0,    0.0,\
  -14000.0, 14000.0,  -20.0]  # -5000

# crystal properties 
a = np.reshape(np.array((1,0,0)),(3,1))  # primitive vector "a"
b = np.reshape(np.array((np.cos(60*np.pi/180),np.sin(60*np.pi/180),0)),(3,1))  # primitive vector "b" for p6
c0 = np.reshape(np.array((0.5,0.5*np.tan(30*np.pi/180),0)),(3,1)) # vector pointing from the UC corner to the oen of the C3 centers in a p6 symmetry
c1 = 2*c0 # second C3 center in a p6 (also rotated by 180 around Z axis)
b0 = a/2 # 90 rotation
b1 = b/2 # 150 rotation
b2 = (a+b)/2 # 30 rotation 
UC_spacing = int(537*2*np.cos(30*np.pi/180)/3)


# colors selection 
compAcolor = [218, 112, 214]  ; cmd.set_color( "colorCompA",compAcolor)
compBcolor = [127,255,212]    ; cmd.set_color( "colorCompB",compBcolor)
compAcolorReleased = [158, 52, 154]  ; cmd.set_color( "colorCompAr",compAcolorReleased)
compBcolorReleased = [77,205,152]    ; cmd.set_color( "colorCompBr",compBcolorReleased)




#### begin simulations
### simulation intervals: determining new positions and translating to those positions 
N=10 # system size (number of A component in the system)
FPS = 30 # frames per second
prob=0.1	# ditachment probability
scale = 15  # movment distance per simulation cycle
XYZrot = np.array([60,0,0])  # rotation target during simulation 
XYZrot2 = np.array([0,720,0]) # rotation target during simulation 
simLen = 18# length of resulting simulation in sec [sec ]      


# build system 
objListAll = setUpAssembly(N,UCtype="hexagon32",colors=["colorCompA","colorCompB"])
colorStart = [0.0,0.0,0.0] ; colorCurrent = colorStart
cmd.set_color( "colorBG",colorStart ) # staring background color
colorEnd = list(np.array([50,30,8])/255)  # dark borwn
current_view = view0

# dump first image
count = 0 
print("p2print: ",os.path.join(p2f,str(count).zfill(4)+".png") )
cmd.set_view(current_view)  # set starting position 
cmd.bg_color("colorBG")
cmd.png(os.path.join(p2f,str(count).zfill(4)+".png"), width=pixelsX, height=pixelsY,ray=1)




for n in range(int(simLen*FPS)):
	#print("n: ",n)
	if n > 1*FPS:
		objListAll = simulationCycle(objListAll,prob,scale,a,b,UC_spacing)
	
	current_view = movCamera(view0,current_view,viewFinal,n,0.5,4,FPS)  # (view,viewFinal,current_itter,start_time,time,FPS)
	current_view = rotCamera(current_view,XYZrot,n,1,3.5,FPS) #rotCamera(view,XYZrot,current_itter,start_time,time,FPS)
	current_view = rotCamera(current_view,XYZrot2,n,3.5,simLen-3.5,FPS)
	colorCurrent = changeColor(colorStart,colorEnd,n,0.5,4,FPS)
	cmd.set_color( "colorBG",colorCurrent ) ; cmd.bg_color("colorBG")
	view0_np = np.array(view0) ; view0_np[:9] = np.array(current_view)[:9]
	view0 = list(view0_np)

	cmd.set_view(current_view)
	cmd.png(os.path.join(p2f,str(n+1).zfill(4)+".png"), width=pixelsX, height=pixelsY ,ray=1)


























