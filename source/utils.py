

# import dependencies
import os, sys, copy, math, itertools, time
import numpy as np
import scipy as sc
from scipy.spatial.transform import Rotation as R
from scipy.spatial import distance
import pymol
from pymol import cmd
from pymol.cgo import *
from pymol.cgo import COLOR, SPHERE, CYLINDER, BEGIN, TRIANGLE_STRIP, NORMAL, VERTEX, END, ALPHA
from pymol.importing import filename_to_format
from Bio.PDB.PDBParser import PDBParser



class compPymolObj:
    name = "" # Name for the computational part
    objName = "" # Name for a pymol object to interact with 
    compType = "" # component type [A/B]
    Sym = 0  # order of the rotation symmetry 
    loc = np.empty((3,1))  # COM of onject
    orientation = np.empty((3,3)) # orientation of the object 
    bound = True  # state of object, is part of an array 
    locked = True # state of object within the array, True if all directions are occupied by other components
    color = [1.0,1.0,1.0]
    rotate = np.zeros((3,1)) # direction of rotating [x,y,z], default is 0
    trajectory = np.zeros((3,1)) # direction of movment [x,y,z], default is 0

    def __init__(self, name,objName,compType,sym,loc,Zangle ):
        self.name = name
        self.objName = objName
        self.compType = compType
        self.sym = sym
        self.loc = np.reshape(loc[:3],(3,-1))
        self.orientation = R.from_euler('xyz', [0,0,Zangle], degrees=True) # to read: np.round(rr.as_matrix(),4)
        
    
    def get_color(self):
        if self.compType=="B":
            return [n/255 for n in [127,255,212]] # Aquamarine 
        if self.compType=="A":
            return [n/255 for n in [218, 112, 214]] # orchid
        else:
            return self.color



def setUpAssembly(N,UCtype="hexagon32",colors=["colorCompA","colorCompB"]):
	if UCtype=="hexagon32":
		# unit cell promitive vectors and geometical parameters 
		a = np.reshape(np.array((1,0,0)),(3,1))  # primitive vector "a"
		b = np.reshape(np.array((np.cos(60*np.pi/180),np.sin(60*np.pi/180),0)),(3,1))  # primitive vector "b" for p6
		c0 = np.reshape(np.array((0.5,0.5*np.tan(30*np.pi/180),0)),(3,1)) # vector pointing from the UC corner to the oen of the C3 centers in a p6 symmetry
		c1 = 2*c0 # second C3 center in a p6 (also rotated by 180 around Z axis)
		b0 = a/2 # 90 rotation
		b1 = b/2 # 150 rotation
		b2 = (a+b)/2 # 30 rotation 
		UC_spacing = int(537*2*np.cos(30*np.pi/180)/3)

		### generate lattice parameters array for the A component (C3)
		coefArray = genA_RotCentersInd(N,a,b,c0,UC_spacing,shape='round')  # note that coefArray[:,:,x] also determine the rotation of the ASU
		Adata = genA_RotCentersPos(coefArray,a,b,c0,UC_spacing)[1] # this was use to generate the lattice fro ma single ASU object
		
		### generate lattice parameters array for the B component (C2)
		coefArrayB = genB_RotCentersInd(coefArray)
		Bdata = genB_RotCentersPos(coefArrayB,a,b) # here distances are still normalized (need to multiply with UC_spacing)
		Bdata[:,:2] = Bdata[:,:2]*UC_spacing
		
		#### databases to a list of objects 
		objListA = [eval("compPymolObj")("A"+str(n[0]).zfill(4),"A"+str(n[0]).zfill(4),"A",3,n[1][:3],n[1][4])  for n in enumerate(Adata)]
		objListB = [eval("compPymolObj")("B"+str(n[0]).zfill(4),"B"+str(n[0]).zfill(4),"B",2,n[1][:3],n[1][4])  for n in enumerate(Bdata)]
		objListAll = objListA+objListB
		cmd.hide(representation="everything",  selection="Acomp")
		cmd.hide(representation="everything",  selection="Bcomp")
		### initiation - generating and positioning all the objects 
		for objComp in objListAll:  # objListAll: list of objects of class compPymolObj instances to allow visualization  
			#print("ASU location: ",latticePos[k])
			createNewObj(objComp) # create a new pymol object namde (objComp.objName) from template named objComp.compType 
			
			#cmd.color("blue",  objNew+" and chain A+C+E+G+I+K")
			#cmd.show( representation="dots",  selection=objNew )
			cmd.translate(list(objComp.loc),  selection = objComp.objName,camera = 0)
			rotationAngleAroundZ = np.round(objComp.orientation.as_euler('xyz')[2]*180/np.pi,4)
			cmd.rotate([0,0,1],rotationAngleAroundZ,selection=objComp.objName ,camera=0,origin=list(objComp.loc))
			cmd.color('colorCompA',objComp.objName) if objComp.compType=="A" else cmd.color('colorCompB',objComp.objName)
		return objListAll



def simulationCycle(objList,prob,scale,a,b,UC_spacing):
	objList = isLock(objList,a,b,UC_spacing)
	objList = isBound(objList,prob,a,b,UC_spacing)
	objList = setTrajectory(objList,a,b,UC_spacing) if np.ceil(np.random.rand(1)-0.7) else objList
	objList = setOrientation(objList)
	for n,objComp in enumerate(objList):
		#if objComp.locked==False:
		#	cmd.color('yellow',objComp.objName)
		if objComp.bound==False:
			if objComp.compType=="A":
				cmd.color('colorCompAr',objComp.objName)
			if objComp.compType=="B":
				cmd.color('colorCompBr',objComp.objName)
	
		objComp.loc = objComp.loc+objComp.trajectory*scale
		cmd.translate(list(objComp.trajectory*scale),  selection = objComp.objName,camera = 0)
	return objList


	

def createNewObj(obj):
	'''
	obj is instance of class compPymolObj
	Object template name is stored in (obj.compType) 
	object name (obj.objName) is the name allowing to interact with pymol 
	'''
	print("Creating obj: ", obj.objName )
	cmd.create(obj.objName,obj.compType+"comp")
	cmd.hide(representation="everything",  selection=obj.objName)
	cmd.show(representation="ribbon",  selection=obj.objName)
	#cmd.set_color( "X%s"%(str(ChainNum)),colorList[ChainNum] ) # p2pdb[0][:-4] and "chain " c
	#cmd.color( "X%s"%(str(ChainNum)), "%s and chain %s"%(ob,c) ) # p2pdb[0][:-4] and "chain " c


def genA_RotCentersInd(N,A,B,C0,UC_spacing,shape='round'):  # round
    '''
    Generating the list of lattice positions
    protocol:
    1. Given a number of latice positions devided by two (to ASUs per unit cell) find the upper bound number of unit cells needed to generae the lattice
    2. Generate list of dtype'i,i,i'
    3. select the number of unit cells needed of the provided list
    
    Inputs: N - ASUs #
            A,B - normalized lattice primitive vectors np.array(3,1)
            C0 - normalized vector to the C3 lower left ASU center from the lower left corner. 2*C0 would be to the second ASU
            shape - shape of the lattice, accepting string of either ['round','Hexagonal']
            
    return: ASUsIndices - np.array(N,4) - [N (number of points), aIndex[int],bIndex[int],c0Index[0,1,2], dist ] 
            where a/bIndeces are the indices of the primitive vector and c0 is the index number of the rotation center,
            and dist is the distance of each point from the center of axes.
    '''
    N_UC = np.ceil(N/2).astype(int) # number of unit cells
    latticeShape = ['round','Hexagonal']
    if shape not in latticeShape:           # n is the max calue of the primitive cell indices 
        print("Input lattice shape is not specified correctly, options are ['round','Hexagonal'], returning a round lattice shape")
        shape = 'round'
    if shape=='round':
        n = np.ceil(np.power(N_UC,0.5)*2/np.sqrt(np.pi)+5).astype(int)   
    if shape=='Hexagonal':
        n = np.ceil(np.power(N_UC,0.5)).astype(int)
        
    centerIndicesVec = np.arange(n,dtype=int)- int((n-1)/2)   # generate a centered vector for the UC indices 
    UCindicesPairs = itertools.product(centerIndicesVec, repeat=2)                      # geenrate pairs of indices to specify each unit cell
    UC_ASUsIndices = np.array([(np.append(np.array(n),0), np.append(np.array(n),1)) for n in  UCindicesPairs])  # generate a double array, one for each ASU
    UC_ASUsIndices = np.reshape(UC_ASUsIndices,[UC_ASUsIndices.shape[0]*2,1,3])  # reshaping to obtain the arrays of locations 
    UC_ASUsIndices = UC_ASUsIndices.flatten()
    UC_ASUsIndices = np.reshape(UC_ASUsIndices,(int(UC_ASUsIndices.shape[0]/3),3))
    
    # ASU selection to create lattice in a given shape 
    if shape not in latticeShape:           # n is the max calue of the primitive cell indices 
        print("Input lattice shape is not specified correctly, returning a round lattice shape")
        shape = 'round'
    if shape=='round':
        p = genA_RotCentersPos(UC_ASUsIndices,A,B,C0,UC_spacing)[0]  # generate physical location (not index) of each position 
        d = np.reshape(np.sum(np.power(p,2),axis=1),[p.shape[0],1])  # distance from center of axis of each point  
        UC_ASUsIndices = np.concatenate((UC_ASUsIndices,d),axis=1)  # append distant as a forth column  
        ind=np.argsort(UC_ASUsIndices[:,-1]) # sort according to distance from the center 
        return UC_ASUsIndices[ind][:N,:]   # UC_ASUsIndicesSelected
    if shape=='Hexagonal':
        UC_ASUsIndices = np.concatenate((UC_ASUsIndices,np.reshape(np.sum(UC_ASUsIndices,axis=1),[UC_ASUsIndices.shape[0],1] )),axis=1)  # sum all the column into a forth column
        ind=np.argsort(UC_ASUsIndices[:,-1])
        UC_ASUsIndicesSelected=UC_ASUsIndices[ind][:N,:]
        
    return UC_ASUsIndicesSelected


def genA_RotCentersPos(arrI,a,b,c0,UC_spacing):
    '''
    input: np.array(N,4) [N (number of points), [aIndex[int],bIndex[int],c0Index[0,1],dist] ]
           a - primitive vector a 
           b - primitive vector b
           c0 - vector from the unit cell origin to the first x1 and second x2 position of the ASU object
                (this is correct for the C3 centers in a P6 plane symmetry only)
    returns: np.array(N,3) N-number of points, (x,y,z) normalized 
    '''
    cImage = np.zeros((3,1)) # constraucting 3D matrix of primitive vectors
    primitiveVecMatrix = np.concatenate((a,b,cImage),axis=1).T
    arr3D = np.concatenate((arrI[:,:2],np.zeros((arrI.shape[0],1))),axis=1) # 3D indeces matrix
    normPosition = np.matmul(arr3D,primitiveVecMatrix)+c0.T
    
    a1 = np.where(arrI[:,2]==1)[0] ; normPosition[a1]=normPosition[a1]+c0.T
    AposType2D = np.reshape(arrI[:,2],(arrI[:,2].shape[0],1)) # creating a 2D array of the column with the rotation symmetry info
    Adata = np.concatenate((normPosition*UC_spacing,AposType2D,AposType2D),axis=1 ) # each line data is [x,y,z,centerType] next add rotation over Z]
    a0 = np.where(arrI[:,2]==0)[0] ; a1 = np.where(arrI[:,2]==1)[0]
    Adata[a0,4]=30 ; Adata[a1,4]=90
    return normPosition,Adata


def genB_RotCentersInd(indArray):
    '''
    Generates the indices matrix of the second components (c2) in a p6 plane
    input: np.array(N,4) - [N (number of points), [aIndex[int],bIndex[int],c0Index[0,1],dist] ]    
    returns: np.array(N,3) - [N (number of points), aIndex[int],bIndex[int],c0Index[0,1,2] ] 
        where a/bIndeces are the indices of the primitive vector and c0 is the index number of the rotation center
    To create a single object in each position we selected to create two object for each previous object with index==0
    and an additional single object to a previous object with index==1
    '''
    a1 = np.where(indArray[:,2] == 0)[0] ; a2 = np.where(indArray[:,2] == 1)[0]
    b0 = indArray[a1,:3] ; b1 = b0.copy() ; b1[:,2]=1
    b2 = indArray[a2,:3] ; b2[:,2]=2
    return np.concatenate((b0,b1,b2), axis=0)
    

def genB_RotCentersPos(arr,a,b):
    '''
    input: np.array(N,3) [N (number of points), [aIndex[int],bIndex[int],c0Index[0,1]] ]
           a - primitive vector a 
           b - primitive vector b
    returns: np.array(N,5) N-number of points, (x,y,z, type, Zrotation) normalized 
    '''
    arr3D = np.concatenate((arr[:,:2],np.zeros((arr.shape[0],1))),axis=1) # 3D mindex matrix
    cImage = np.zeros((3,1)) # constraucting 3D matrix of primitive vectors
    primitiveVecMatrix = np.concatenate((a,b,cImage),axis=1).T
    latticePos = np.matmul(arr3D,primitiveVecMatrix) # UC origin position, final position is determined by the type number
    
    a0 = np.where(arr[:,2]==0)[0]  # all locations where element index ==0
    latticePos[a0] = latticePos[a0]+a.T/2 # appending the local vector
    a1 = np.where(arr[:,2]==1)[0] ; latticePos[a1] = latticePos[a1]+b.T/2
    a2 = np.where(arr[:,2]==2)[0] ; latticePos[a2] = latticePos[a2]+(a+b).T/2
    dataB = np.concatenate((latticePos,np.zeros((latticePos.shape[0],2))),axis=1)
    dataB[:,3]=arr[:,2]  # adding center type on fourth column
    dataB[a0,4]=90 ; dataB[a1,4]=150 ; dataB[a2,4]=30
    return dataB


def getUCp0(N,a,b,S):
    '''
    N - numer of ASUs to be imaged 
    a,b - primitive vectors
    S - unit cell spacing
    '''
    coefArray = genA_RotCentersInd(N)
    UCp0List = np.unique(coefArray[:,:2],axis=0)
    UCp0Locs = np.dot(UCp0List,[a,b])*S  
    UCp0Locs = np.append(UCp0Locs, np.zeros((UCp0Locs.shape[0],1)), axis=1)
    return UCp0List,UCp0Locs


def objList2npXYZ(objList):
    '''
    Function recived list of instances of class compPymolObj
    returns: np.array(Nx3) where N is the length of objList (number of instances)
    '''
    return np.reshape(np.array([n[1].loc for n in enumerate(objList)]),(len(objList),3) )


def getAdjacentLocations(obj,a,b,UC_spacing):
    '''
    Function recived list of instances of class compPymolObj, primitive vectors and spacing scalar 
    returns: np.array(Nx3) where N is the number of direction an object connect, and 3 fro XYZ 
    '''
    x0 = obj.loc # start point of vector
    sym = obj.sym # symmetry of object around the Z axis 
    rotationVec = np.arange(0,360,360/sym) # binding direction assuming it is aligned with the z and x axis 
    otient0 = obj.orientation # rotation object
    # creating list of rotation object for each binding direction (number depends on symmetry order)
    rotationObjList = [otient0*R.from_euler('xyz', [0,0,Zrot], degrees=True) for Zrot in rotationVec]
    A2Bdist = np.linalg.norm(a)/2*np.tan(30*np.pi/180)*UC_spacing # distant between COM of components in assembly state
    vec = np.array((A2Bdist,0,0)) # initializing 1D vector
    return np.reshape(np.array([x0+np.reshape(rotObj.apply(vec),(3,1)) for rotObj in rotationObjList]),(len(rotationObjList),3))
    

def boundInterfaceCount(obj,objList,a,b,UC_spacing):
    '''
    task: return number of interfaces that are in contact
    input: obj to examine
           list of all objects in the system
    output: number of contacts 
            lock state, true if number of contact is equal to sym order
    '''
    V = getAdjacentLocations(obj,a,b,UC_spacing) # get the object adjucent locations in its current state
    db = objList2npXYZ(objList) # get all the system current locations 
    A2Bdist = np.linalg.norm(a)/2*np.tan(30*np.pi/180)*UC_spacing # distant between COM of components in assembly state
    cutOffRadius = A2Bdist/5 # cut of distance for COM position is arbitrarily set to 1/5 of the system distance
    minDist = np.min(distance.cdist(V, db),axis=1)
    
    numOfContacts = sum([1  for n in np.min(distance.cdist(V, db),axis=1) if n<cutOffRadius])
    sym = obj.sym
    lock=1 if numOfContacts==sym  else 0 
    return numOfContacts,lock
        
    
def unBlockedTrajectory(obj,objList,a,b,UC_spacing):
    '''
    task: return number of interfaces that are in contact
    input: obj to examine
           list of all objects in the system
    output: number of contacts 
            lock state, true if number of contact is equal to sym order
    '''
    V = getAdjacentLocations(obj,a,b,UC_spacing) # get the object adjucent locations in its current state
    db = objList2npXYZ(objList) # get all the system current locations 
    A2Bdist = np.linalg.norm(a)/2*np.tan(30*np.pi/180)*UC_spacing # distant between COM of components in assembly state
    cutOffRadius = A2Bdist/5 # cut of distance for COM position is arbitrarily set to 1/5 of the system distance
    minDist = np.min(distance.cdist(V, db),axis=1)
    vdirect = V[np.where(minDist>cutOffRadius)[0]]
    #print(vdirect,vdirect.shape,vdirect[0])
    return np.reshape(np.round(vdirect[0,:]/np.linalg.norm(vdirect[0,:]),4),(3,1))
    
    
def isLock(objList,a,b,UC_spacing):
	'''
	Updates the class instances with the system relevant locked variable
	input: all the intances of class compPymolObj
	ouput: updated .locked value (default was true)
	A component is considered locked if all its interfaces are in contact, here the calculation means
	That in each potantial position the system can find another component COM in 2nm from the expected point
	'''
	for obj in objList:
		n,l = boundInterfaceCount(obj,objList,a,b,UC_spacing)   
		obj.locked=True if l else False
	return objList
 

def isBound(objList,prob,a,b,UC_spacing):
    '''
    Updates the class instances with the system relevant bound variable. 
    A component is considered bound if at list one of its interfaces is in contact. 
    input: all the intances of class compPymolObj
    ouput: updated .bound value (default was true)
    '''
    for obj in objList:
        if obj.locked==False:
            if obj.bound==True: 
                n,l = boundInterfaceCount(obj,objList,a,b,UC_spacing)
                P = np.power(prob,n)
                boundState = np.random.choice((False,True), p=[P, 1-P])
                obj.bound=boundState 
        if obj.locked==True:
            obj.bound=True
    return objList
       

def setTrajectory(objList,a,b,UC_spacing):
    '''
    Updates the class instances with trajectory value pending the bound/lock values and system COM values
    A component can have a non zero trajectory value only if it is not locked. 
    input: objList all the system intances of class compPymolObj, prob [0,1] protbability for a bound unlocked component to move
    ouput: updated .tragectory value (default was [0,0,0])
    '''
    for N,obj in enumerate(objList):
         if obj.bound==False:
            n,l = boundInterfaceCount(obj,objList,a,b,UC_spacing)
            if not n: 
                vec  = np.random.rand(3)-0.5 
                obj.trajectory = np.reshape(vec/np.linalg.norm(vec),(3,1)) #  random move
            else: # case obj.bound==False but is still in the vicinity of components 
                obj.trajectory = unBlockedTrajectory(obj,objList,a,b,UC_spacing) 
    return objList


def setOrientation(objList):
    '''
    Updates the class instances with orientation matrix - this will occure only for non bound component 
    input: objList all the system intances of class compPymolObj
    ouput: updated .orientation value with randomely generated orientation matrix 
    
    '''
    for N,obj in enumerate(objList):
        if obj.bound==False:
            anglesVector = (np.random.rand(3)-0.5)*30 #  degrees 
            RotObj = R.from_euler('xyz',anglesVector,degrees=True) # arbitrary rotation object around each axis in the range of [-45,45] degrees
            obj.orientation = obj.orientation*RotObj 
            cmd.rotate([1,0,0],anglesVector[0],selection=obj.objName ,camera=0,origin=list(obj.loc))
            cmd.rotate([0,1,0],anglesVector[1],selection=obj.objName ,camera=0,origin=list(obj.loc))
            cmd.rotate([0,0,1],anglesVector[2],selection=obj.objName ,camera=0,origin=list(obj.loc))
    return objList


def cgoCylinder3(p0=[0,0,0],v0=[1,0,0],inputAng=0, L=20, radius=0.5 ,color=[1,0,0]  ,name=''):
	p0 = np.array(p0)
	V0 = np.array(v0)
	r = R.from_euler('z', inputAng, degrees=True)
	p1 = np.around(R.apply(r,V0), decimals=6)*L+p0
	PP0 = list(p0) ; PP1=list(p1)
	obj  = [cgo.CYLINDER] + PP0 + PP1 + [radius] + color + color
	if not name:
		name = cmd.get_unused_name('cylinderObj')
	cmd.load_cgo(obj, name)
	return name,obj, p1


def drowUnitCell3(Type="Hexagonal", UCp0=[0,0,0],UCa=150,UCb=150,V0=[1,0,0] ,UCaAng=0,UCbAng=45,UCc=[218/256,165/256,32/256], UCr=1):
	'''
	Type="Hexagonal", UCp0=UCp0,UCa=UCa,UCb=UCb,UCaAng=UCaAng,UCbAng=UCbAng,UCc=UCc,UCr=UCr
	'''
	UCp0 = np.array(UCp0)
	V0 = np.array(V0)
	UCtypes = ["Hexagonal","Oblique","Rectangular","square","Rhombic"]
	#V0 = [1,0,0] # Default orientation (facing the z direction) for plane in XY plane
	if Type not in UCtypes:
		print(Type ,"unit cell does not exist") ; sys.exit()
	elif Type=="Hexagonal":
		UCbAng = 60  ; UCb=UCa  # input angles
	elif Type=="square":
		UCbAng = 90  ;  UCb=UCa  # input angles
	elif Type=="Oblique":
		None   # input angles
	elif Type=="Rectangular":
		UCbAng = 90    # input angles
	elif Type=="Rhombic":
		UCb=UCa  # input angles

	name1,obj1, p1 = cgoCylinder3(p0=UCp0,v0=V0,inputAng=UCaAng,L=UCa,radius=UCr,color= UCc)
	name2,obj2, p2 = cgoCylinder3(p0=UCp0,v0=V0,inputAng=UCbAng,L=UCb,radius=UCr,color= UCc)
	name3,obj3, p3 = cgoCylinder3(p0=p2,v0=V0,inputAng=UCaAng,L=UCa,radius=UCr,color= UCc)
	name4,obj4, p4 = cgoCylinder3(p0=p1,v0=V0,inputAng=UCbAng,L=UCb,radius=UCr,color= UCc)
	objList = [name1,name2,name3,name4]
	print('optional unnit cells: ["Hexagonal","Oblique","Rectangular","square","Rhombic"]' )
	return UCp0,p1,p2,p3,objList


def drow2Dlattice(gName="",Type="Hexagonal",UCp0=[0,0,0],UCa=310,UCb=100,UCaAng=0,UCbAng=60,UCc=[218/256,165/256,32/256],UCr=1,SE=1,SEl=15,SEr=1.2):
	'''
	Start input some documentation here
	This fucnction generate a pymol set of CGOs to illustrate the p632 plane symmetry 
	inputs: 
	UCp0 - origin of unit cell list/array of dimension 3 
	UCa - first unit cell primitive vector length (unit cell spacing) in Agstrom
	UCb - second unit cell primitive vector length (unit cell spacing) in Agstrom, for some unit cell types this value is predetermined as equal to the first primitive vector
	UCaAng - direction of the first primitive vector , by default along the x axis
	UCbAng - angle of the second primitive vector in x,y plane vs. the first primitive vector, this is a variables only for unit cells of type Oblique and Rhombic
	UCc - cgo cylinders color (currently type of orange)
	UCr - cgo cylinders radio=us. 
	SE = boolian, drow symmetry elements 
	SEl - diameter of rotation symmetry object to be printed 
	SEr - cgo cylinders radius belonging to the rotation symmetry objects

	Output:
	gName = name of group all the objects are clustered under  
	objList - list of objects geerated  (this allows pther functions to modify or delete produced selections. )
	np.array(UCp0),p1,p2,p3 - Unitcell corners position in x,y,z list format
	
	'''
	SEc1 = [1,0,0] ; SEc2 = [0,1,0] ; SEc3 = [0,1,1] ; SEc4 = [0,0,1]  # colors of symmetry operations 
	UCp0 = np.array(UCp0)
	V0 = np.array([1,0,0]) # Default orientation (facing the x direction) in the plane in XY plane 
	UCp0,p1,p2,p3,objList1 = drowUnitCell3(Type="Hexagonal", UCp0=UCp0,UCa=UCa,UCb=UCb,V0=V0,UCaAng=UCaAng,UCbAng=UCbAng,UCc=UCc,UCr=UCr)
	if SE:  # default to drow the sym elements 
		nfold = 6 
		objList11=drowRect(UCp0,SEl,nfold,SEr,SEc1) ; objList12=drowRect(p1,SEl,nfold,SEr,SEc1) ;objList13=drowRect(p2,SEl,nfold,SEr,SEc1) ;objList14=drowRect(p3,SEl,nfold,SEr,SEc1)		
		nfold = 3 ; P3_1 = 2/3*np.array(UCp0)+1/3*p3 ; P3_2 = 1/3*np.array(UCp0)+ 2/3*p3
		objList21=drowRect(P3_1,SEl,nfold,SEr,SEc2) ; objList22=drowRect(P3_2,SEl,nfold,SEr,SEc2)
		nfold = 2 ; SEr=3;SEl=10; P2_1 = (np.array(UCp0)+p1)/2 ; P2_2 = (np.array(UCp0)+p2)/2 ; P2_3 = (p1+p3)/2 ;P2_4 = (p2+p3)/2 ; P2_5 = (np.array(UCp0)+p3)/2 
		objList31=drowRect(P2_1,SEl,nfold,SEr,SEc4) ; objList32=drowRect(P2_2,SEl,nfold,SEr,SEc4) ; objList33=drowRect(P2_3,SEl,nfold,SEr,SEc4); objList34=drowRect(P2_4,SEl,nfold,SEr,SEc4) ; objList35=drowRect(P2_5,SEl,nfold,SEr,SEc4)
	
	### preparing pumol CGO selection output groups 
	objList =  objList11+objList12+objList13+objList14+objList21+objList22+objList31+objList32+objList33+objList34+objList35
	### list of unique names not in use in the pymol session
	gName = cmd.get_unused_name('p6')
	frameName = cmd.get_unused_name('frame')
	elementsName = cmd.get_unused_name('elements')
	### creating two sub groups of the unit cell frame and the symmetry elements within
	[ cmd.group(frameName, members=name , action='add', quiet=1)  for name in objList1]
	[ cmd.group(elementsName, members=name , action='add', quiet=1)  for name in objList]
	### group of the two above
	[ cmd.group(gName, members=name , action='add', quiet=1)  for name in [frameName,elementsName]]
	print("Unit Cell ",gName," Corbers are in: ", np.array(UCp0),p1,p2,p3)
	return np.array(UCp0),p1,p2,p3,objList,gName


def drow2DlatticeFotN_ASUs(N,Type="Hexagonal",UCp0=[0,0,0],UCa=310,UCb=100,UCaAng=0,UCbAng=60,SEl=50,SEr=1.2):
	'''
	Unit cell primitive vectors definition and unit cell spacing
	notes that those definitions void the inputs UCp0=[0,0,0] UCa=310,UCb=100,UCaAng=0,UCbAng=60 something that would need to be considered later in the future
	'''
	a = np.array((1,0))
	b = np.array((np.cos(60*np.pi/180),np.sin(60*np.pi/180)))
	c0 = np.array((0.5,0.5*np.tan(30*np.pi/180)))
	c1 = 2*c0
	UC_spacing = int(537*2*np.cos(30*np.pi/180)/3)

	a0,a1 = getUCp0(N,a,b,c0,UC_spacing)  # a0 is the unit cells id and a1 is the p0 position of each unit cell. 
	print(a0,a1)
	count=0
	for p0 in a1:  # utterating through all ASUs starting positions 
		se=0 if count>0 else 1
		drow2Dlattice(gName="",Type="Hexagonal",UCp0=p0,UCa=310,UCb=100,UCaAng=0,UCbAng=60,UCc=[218/256,165/256,32/256],UCr=1,SE=se,SEl=50,SEr=1.2)
		

### pymol function for create, translate and rotate 
def createAndTranlateObj(obj,count=0,name="ASU"):
	'''
	Take an object and create a new identical object with unique name 
	obg - pymol object
	count - number of the new object for a unique id
	name - blue print for the new objects name
	'''
	objNew = name+str(count).zfill(6)
	print("Creating obj: ", objNew )
	cmd.create(objNew,obj)
	cmd.hide(representation="everything",  selection=objNew)
	cmd.show(representation="ribbon",  selection=objNew)
	return objNew


# getObgSeqLen("ASU")
def getObgSeqLen(obj):
	obj_names = cmd.get_names('objects', 1, obj)
	Chains = cmd.get_chains('%s and (%s)' % (obj_names[0],"all"))
	chainSele = [obj_names[0]+" and chain "+chain for chain in Chains]  # list of selection strings 	
	resList = [[at.resn for at in cmd.get_model(c).atom] for c in chainSele]
	resNum = len(list(itertools.chain(*resList)))
	print("Number of Atoms in the model = ",resNum)
	return resList,resNum


# creatDilutedObj("ASU")
def creatDilutedObj(obj,ratio=0.5, name="ASUbase"):	
	'''
	recieve an object and create an object with diluted atoms base on the ratio input 
	'''
	obj_names = cmd.get_names('objects', 1, obj)
	Chains = cmd.get_chains('%s and (%s)' % (obj_names[0],"all"))
	chainSele = [obj_names[0]+" and chain "+chain for chain in Chains]  # list of selection strings 	
	resList = [[at.resi for at in cmd.get_model(c).atom] for c in chainSele]
	resIndx = list(itertools.chain(*resList))
	
	dilute = np.ceil(1/ratio).astype(int)
	diluteResList = resIndx[0::dilute]
	print('+'.join(diluteResList))
	print("resIndx: " ,len(resIndx))
	print("diluteResList: ",len(diluteResList))

	diluteResSele = '+'.join(diluteResList)
	cmd.create(name,obj_names[0]+" and resi "+diluteResSele)
	print("Creating an object with ",len(diluteResList)," atom from the ",len(resIndx)," atom containing input object"  )
	return name


def drowN_ASUs(N,Type="Hexagonal",shape='round',UCp0=[0,0,0],UCa=310,UCb=100,UCaAng=0,UCbAng=60,name="ASU"):

	# note the 5 lines below are design specific and surpass the default input UCa, UCb, UCaAng, and UCbAng (which can than becaluculated from a,b, and UC-spacing) 
	a = np.array((1,0))
	b = np.array((np.cos(60*np.pi/180),np.sin(60*np.pi/180)))
	c0 = np.array((0.5,0.5*np.tan(30*np.pi/180)))
	c1 = 2*c0
	UC_spacing = int(537*2*np.cos(30*np.pi/180)/3)
	##########################

	#name="ASU"  # for selection purposes 
	coefArray = genA_RotCentersInd(N,a,b,c0,shape='round')  # indices of the UC for each ASU, the third element is the Z orientation (if 1 rotation of 180 degrees is applied)
	latticePos = genA_RotCentersPos(coefArray,a,b,c0)*UC_spacing        #  this is the set of ASUs CM position 

	LofL = [] #  Objects list  
	for n,k in enumerate(range(N)):
		print("ASU location: ",latticePos[k])
		objNew = createAndTranlateObj("ASU",count=k,name="ASU") # create a new object and return name of object 
		cmd.translate(list(latticePos[k]),  selection = objNew,camera = 0)
		if coefArray[k,2]==1:
			cmd.rotate([0,0,1],180,selection=objNew ,camera=0,origin=latticePos[k])
		LofL.append(objNew)
		#LofL.append([s]+LofL[n-1])  if n>0 else LofL.append([s])  # create list of list, each list contain 
	print("Arrays object name: ", LofL)


def get_str(format, selection='(all)', state=-1, ref='',ref_state=-1, multi=-1, quiet=1, _self=cmd):
	'''
	DESCRIPTION
	Like "get_bytes" but return a unicode string.
	'''
	assert format not in ('mmtf',), 'binary format, use get_bytes'
	b = get_bytes(format, selection, state, ref, ref_state, multi, quiet, _self)
	if b is None:
		return None
	return b.decode('utf-8')


def multisave(filename, pattern="all", state=-1,append=0, format='', quiet=1, _self=cmd):
	print(1)
	'''
	DESCRIPTION
	
	"multisave" will save a multi-entry PDB file.
	
	Every object in the given selection (pattern) will have a HEADER and a
	CRYST (if symmetry is defined) record, and is terminated with END.
	Loading such a multi-entry PDB file into PyMOL will load each entry
	as a separate object.
	
	This behavior is different to the "save" command, where a multi-object
	selection is written "flat" to a PDB file, without HEADER or CRYST
	records.
	
	ARGUMENTS
	
	filename = string: file path to be written
	pattern = str: atom selection (before 1.8.4: object name pattern)
	state = int: object state (-1=current, 0=all) {default: -1}
	append = 0/1: append to existing file {default: 0}
	format = str: file format {default: guess from extension, or 'pdb'}
	'''
	_, _, format_guessed, zipped = filename_to_format(filename)
	if zipped:
		raise pymol.CmdException(zipped + ' not supported with multisave')
	if not format:
		format = format_guessed or 'pdb'
	if format == 'pmo':
		raise pymol.CmdException('pmo format not supported anymore')
	if format not in ('pdb', 'cif'):
		raise pymol.CmdException(format + ' format not supported with multisave')
	s = get_str(format, pattern, state, '', -1, 1, quiet, _self)
	if s is None:
		if _self._raising(): raise QuietException
		return DEFAULT_ERROR
	filename = _self.exp_path(filename)
	with open(filename, 'a' if int(append) else 'w') as handle:
		handle.write(s)
	return DEFAULT_SUCCESS


def drowSaveN_ASUs(N,Type="Hexagonal",shape='round',ratio=1,UCp0=[0,0,0],UCa=310,UCb=100,UCaAng=0,UCbAng=60,name="ASU"):

	# note the 5 lines below are design specific and surpass the default input UCa, UCb, UCaAng, and UCbAng (which can than becaluculated from a,b, and UC-spacing) 
	a = np.array((1,0))  # primitive vector "a"
	b = np.array((np.cos(60*np.pi/180),np.sin(60*np.pi/180)))  # primitive vector "b" for p6
	c0 = np.array((0.5,0.5*np.tan(30*np.pi/180))) # vector pointing from the UC corner to the oen of the C3 centers in a p6 symmetry
	c1 = 2*c0 # second C3 center in a p6 (also rotated by 180 around Z axis)
	b0 = a/2 # 90 rotation
	b1 = b/2 # 150 rotation
	b2 = (a+b)/2 # 30 rotation
	UC_spacing = int(537*2*np.cos(30*np.pi/180)/3)
	##########################

	#name="ASU"  # for selection purposes 
	coefArray = genA_RotCentersInd(N,a,b,c0,shape='round')  # indices of the UC for each ASU, the third element is the Z orientation (if 1 rotation of 180 degrees is applied)
	latticePos = genA_RotCentersPos(coefArray,a,b,c0)*UC_spacing        #  this is the set of ASUs CM position 

	if ratio!=1:
		creatDilutedObj(name,ratio,name="ASUbase")
		name = "ASUbase"


	LofL = [] #  Objects list  
	for n,k in enumerate(range(N)):
		print("ASU location: ",latticePos[k])
		objNew = createAndTranlateObj(name,count=k,name=name) # create a new object and return name of object 
		print("objNew: ",objNew)
		cmd.color("blue",  objNew+" and chain A+C+E+G+I+K")
		cmd.show( representation="dots",  selection=objNew )

		cmd.translate(list(latticePos[k]),  selection = objNew,camera = 0)
		if coefArray[k,2]==1:
			cmd.rotate([0,0,1],180,selection=objNew ,camera=0,origin=latticePos[k])
		LofL.append(objNew)
		#LofL.append([s]+LofL[n-1])  if n>0 else LofL.append([s])  # create list of list, each list contain 
	print("Arrays object name: ", LofL)

	### dumping list of selections (each dump will have an additional ASU )
	A = []  # 
	Aname = []
	

def movCamera(view,current_view,viewFinal,current_itter,start_time,time,FPS):
	'''
	Calculate movement of camera 
	input: view - pymol string (18x1) describing the camera properties 
		   viewFinal - destination of camera 
		   current_itter -  simulation  itteration number to get time this is  current_itter/FPS
		   start_time = [sec] at which time during the simulation this should begin
		   time - length of action
		   FPS - number of steps required (time*FPS)
	'''
	current_time = current_itter/FPS
	if(start_time<current_time and current_time<start_time+time):  # this is when function should work
		cameraXYZ = np.array(view[9:12])
		cameraXYZ_Final = np.array(viewFinal[9:12])
		distRatio = (current_time-start_time)/time
		
		currentView = np.array(current_view)
		currentView[9:12]= (cameraXYZ_Final-cameraXYZ)*distRatio+cameraXYZ
		return list(currentView)
	elif current_time<=start_time:
		return view
	else:
		return current_view


def rotCamera(view,XYZrot,current_itter,start_time,time,FPS):
	'''
	Calculate movement of camera 
	input: view - pymol list (18x1) describing the camera properties 
		   XYZrot - np.array(3,1) camera required full rotation around any axes in degrees=True
		   current_itter -  simulation  itteration number to get time this is  current_itter/FPS
		   start_time = [sec] at which time during the simulation this should begin
		   time - length of action
		   FPS - number of steps required (time*FPS)
	'''
	current_time = current_itter/FPS
	if(start_time<current_time and current_time<start_time+time):  # this is when function should work
		cameraOrient = R.from_matrix(np.reshape(np.array(view[:9]),(3,3)))
		currentAngleRotation = R.from_euler('xyz', XYZrot/(time*FPS), degrees=True)
		calcedOrient = np.reshape(np.array((cameraOrient*currentAngleRotation).as_matrix()),(9,))
		view_np = np.array(view)
		view_np[:9] = calcedOrient
		print("turn")
		return list(view_np)
	else:
		return view


def changeColor(colorS,colorE,current_itter,start_time,time,FPS):
	'''
	Calculate movement of camera 
	input: view - pymol list (18x1) describing the camera properties 
		   XYZrot - np.array(3,1) camera required full rotation around any axes in degrees=True
		   current_itter -  simulation  itteration number to get time this is  current_itter/FPS
		   start_time = [sec] at which time during the simulation this should begin
		   time - length of action
		   FPS - number of steps required (time*FPS)
	'''
	current_time = current_itter/FPS
	if(start_time<current_time and current_time<start_time+time):  # this is when function should work
		colorS_np = np.array(colorS)
		colorE_np = np.array(colorE)
		Ratio = (current_time-start_time)/time
		color = list((colorE_np-colorS_np)*Ratio+colorS_np)
		return color
	elif current_time<=start_time: 
		return colorS
	else: 
		return colorE


def drow2Comp_N_ASUs(N,Type="Hexagonal",shape='round',ratio=1,UCp0=[0,0,0],UCa=310,UCb=100,UCaAng=0,UCbAng=60,name="ASU"):

	# note the 5 lines below are design specific and surpass the default input UCa, UCb, UCaAng, and UCbAng (which can than becaluculated from a,b, and UC-spacing) 
	### lattice definitions
	view = [ 1.000000000,    0.000000000,    0.000000000,\
     	 0.000000000,    1.000000000,    0.000000000,\
     	 0.000000000,    0.000000000,    1.000000000,\
     	 0.000000000,    0.000000000, -4313.336425781,\
    	 -0.318004608,    2.596000671,    7.080001354,\
  	 	 4272.680664062, 4353.991699219,  -20.000000000 ]


	view0 = [
	-0.0,    1.0,    0.0,\
	-1.0,   -0.0,   -0.0,\
	-0.0,   -0.0,    1.0,\
	0.0,    0.0, -1350.0,\
	-0.0,    0.0,    0.0,\
	-7000.0, 7000.0,  -20.000000000]

	view_dist = [
	-0.0,    1.0,    0.0,\
	-1.0,   -0.0,   -0.0,\
	-0.0,   -0.0,    1.0,\
	0.0,    0.0, -5000.0,\
	-0.0,    0.0,    0.0,\
	-14000.0, 14000.0,  -20.000000000]

	view_angle = [
	0.5,0.,0.8660254,\
	0. ,  1. , 0. ,\
	-0.8660254,  0.       ,  0.5,\
	0.0,    0.0, -5000.0,\
	-0.0,    0.0,    0.0,\
	-14000.0, 14000.0,  -20.000000000]


	compAcolor = [218, 112, 214]  ; cmd.set_color( "colorCompA",compAcolor)
	compBcolor = [127,255,212]    ; cmd.set_color( "colorCompB",compBcolor)
     



	a = np.reshape(np.array((1,0,0)),(3,1))  # primitive vector "a"
	b = np.reshape(np.array((np.cos(60*np.pi/180),np.sin(60*np.pi/180),0)),(3,1))  # primitive vector "b" for p6
	c0 = np.reshape(np.array((0.5,0.5*np.tan(30*np.pi/180),0)),(3,1)) # vector pointing from the UC corner to the oen of the C3 centers in a p6 symmetry
	c1 = 2*c0 # second C3 center in a p6 (also rotated by 180 around Z axis)
	b0 = a/2 # 90 rotation
	b1 = b/2 # 150 rotation
	b2 = (a+b)/2 # 30 rotation 
	
	UC_spacing = int(537*2*np.cos(30*np.pi/180)/3)
	
	#NN = 50
	### generate lattice parameters array for the A component (C3)
	coefArray = genA_RotCentersInd(N,a,b,c0,UC_spacing,shape='round')  # note that coefArray[:,:,x] also determine the rotation of the ASU
	Adata = genA_RotCentersPos(coefArray,a,b,c0,UC_spacing)[1] # this was use to generate the lattice fro ma single ASU object
	
	### generate lattice parameters array for the B component (C2)
	coefArrayB = genB_RotCentersInd(coefArray)
	Bdata = genB_RotCentersPos(coefArrayB,a,b) # here distances are still normalized (need to multiply with UC_spacing)
	Bdata[:,:2] = Bdata[:,:2]*UC_spacing
	
	#### databases to a list of objects 
	objListA = [eval("compPymolObj")("A"+str(n[0]).zfill(4),"A"+str(n[0]).zfill(4),"A",3,n[1][:3],n[1][4])  for n in enumerate(Adata)]
	objListB = [eval("compPymolObj")("B"+str(n[0]).zfill(4),"B"+str(n[0]).zfill(4),"B",2,n[1][:3],n[1][4])  for n in enumerate(Bdata)]
	objListAll = objListA+objListB

	cmd.hide(representation="everything",  selection="Acomp")
	cmd.hide(representation="everything",  selection="Bcomp")

	### initiation - generating and positioning all the objects 
	for objComp in objListAll:  # objListAll: list of objects of class compPymolObj instances to allow visualization  
		#print("ASU location: ",latticePos[k])

		createNewObj(objComp) # create a new pymol object namde (objComp.objName) from template named objComp.compType 
		

		#cmd.show( representation="dots",  selection=objNew )
		cmd.translate(list(objComp.loc),  selection = objComp.objName,camera = 0)
		rotationAngleAroundZ = np.round(objComp.orientation.as_euler('xyz')[2]*180/np.pi,4)
		cmd.rotate([0,0,1],rotationAngleAroundZ,selection=objComp.objName ,camera=0,origin=list(objComp.loc))
		cmd.set_view(view)

		cmd.color('colorCompA',objComp.objName) if objComp.compType=="A" else cmd.color('colorCompB',objComp.objName)

	### simulation intervals: determining new positions and translating to those positions 
	prob=0.5	
	scale = 15
	objListAll = isLock(objListAll,a,b,UC_spacing)
	objListAll = isBound(objListAll,prob,a,b,UC_spacing)
	objListAll = setTrajectory(objListAll,a,b,UC_spacing)
	for objComp in objListAll:
		#if objComp.locked==False:
		#	cmd.color('yellow',objComp.objName)
		#if objComp.bound==False:
		#	cmd.color('red',objComp.objName)
		objComp.loc = objComp.loc+objComp.trajectory*scale
		cmd.translate(list(objComp.trajectory*scale),  selection = objComp.objName,camera = 0)

	#for k in range(10):
	#	objListAll = isLock(objListAll,a,b,UC_spacing)
	#	objListAll = isBound(objListAll,a,b,UC_spacing)
	#	objListAll = setTrajectory(objListAll,prob,a,b,UC_spacing)
	#	for obj in objListAll:
	#		obj.loc = obj.loc+obj.trajectory*15    
	#		#positions = objList2npXYZ(objListAll)
	#	for obj in objListAll:
	#		cmd.translate(list(obj.loc),  selection = obj.objName,camera = 0)
	#	time.sleep(5)



def saveModel(sele,name ):

	print(1)
	# use: multisave ASU6.pdb, AAA






####################################################################################################################################################
##################################### utils for video preparation #####################################
####################################################################################################################################################
import cv2
import progressbar as pb
from PIL import Image, ImageDraw, ImageFont
from time import strftime
from time import gmtime
import textwrap




def text2video(frames,text,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS=30,t=[0,0]):
    img = Image.new('RGB', (imageSize[0], imageSize[1] ), color = bg_color)
    fnt = ImageFont.truetype(fontSource, fontSize)
    d = ImageDraw.Draw(img)
    d.multiline_text((5, 10), text, fill=textColor, font=fnt)
    #d.text((3,10), text, font=fnt, fill=textColor)
    imgNP = np.array(img)
    
    ### dimensions 
    pixelsY = frames[0].shape[0]
    pixelsX = frames[0].shape[1]
    yStart = position[0]
    yEnd = yStart+imageSize[1]
    xStart = position[1]
    xEnd = xStart+imageSize[0]

    ### timing 
    for k in range(len(frames)):
        if (t[0]==0 and t[1]==0):
            frames[k][yStart:yEnd,xStart:xEnd,:] =  np.clip(frames[k][yStart:yEnd,xStart:xEnd,:].astype('int')+imgNP.astype('int'),0,255).astype('uint8')
        elif (k>FPS*t[0] and k<FPS*t[1]):
            frames[k][yStart:yEnd,xStart:xEnd,:] =  np.clip(frames[k][yStart:yEnd,xStart:xEnd,:].astype('int')+imgNP.astype('int'),0,255).astype('uint8')
    return frames


def clock2video(frames,textColor,bg_color,imageSize,fontSize,fontSource,position,FPS=30):
    TotalTime = strftime("%M:%S", gmtime(len(frames)/FPS))
    ### dimensions 
    pixelsY = frames[0].shape[0]
    pixelsX = frames[0].shape[1]
    yStart = position[0]
    yEnd = yStart+imageSize[1]
    xStart = position[1]
    xEnd = xStart+imageSize[0]

    ### timing 
    for k in range(len(frames)):
        CurrentTime = strftime("%M:%S", gmtime(k/FPS))
        img = Image.new('RGB', (imageSize[0], imageSize[1] ), color = bg_color)
        fnt = ImageFont.truetype(fontSource, fontSize)
        d = ImageDraw.Draw(img)
        d.rectangle([0,0,imageSize[0]-1,imageSize[1]-1], fill=None, outline=(255,255,255))
        d.text((6,5), CurrentTime+" of "+TotalTime, font=fnt, fill=textColor)
        dNP = np.asarray(img).copy().astype(np.uint64)
        timeRatio = int(k/len(frames)*imageSize[0])
        dNP[1:-1,1:timeRatio,:] += 255
        np.place(dNP, dNP>255,0)
        dNP = dNP.astype(np.uint8)
        #new_im = Image.fromarray(dNP)
        #imgNP = np.array(img)
        frames[k][yStart:yEnd,xStart:xEnd,:] =  dNP
    return frames


















