#IF YOU USE THIS SOFTWARE IN YOUR ACADEMIC RESEARCH PLEASE CITE THE PAPER WHERE IT WAS PRESENTED:
#Pietro Speroni di Fenizio: "The Lattice of Chemical Organisations" in Proceedings of European Conference of Artificial Life (ECAL) 2015. 20-24 July, 2015. York. In print.

#in order of appearance :-)
import sys
import shelve
import math
import copy
import itertools


LogOn=False
ShowIntersection=False#True #False
ShowUnion=False#True #False



###########################################
######## General Uses #########################
###########################################


def sortpairs(a,b):
	"""Sorts two list/... based upon their first element """
	if a[1]<b[1]:		return -1
	elif a[1]==b[1]:	return 0
	return 1

class SymetricTable(object):
	"a class to create a table that is symetric. Given a tuple it will order it before checking"
	table={}
	def __init__(self): 
		self.table={}
	def orderTuple(self,x):
		y=list(x)
		y.sort()
		return tuple(y)
	def __contains__(self,key):
		orderedkey=self.orderTuple(key)
		return (orderedkey in self.table)
	def __getitem__(self,key):
		orderedkey=self.orderTuple(key)
		return self.table[orderedkey]		
	def get(self,key,ReturnNotThere=None):
		orderedkey=self.orderTuple(key)
		try:
			return self.table[orderedkey]
		except KeyError:
			return ReturnNotThere
	def __setitem__(self,key,value):
		self.table[self.orderTuple(key)]=value
		return
	def __len__(self): 
		return len(self.table)
	


def LogThis(texttolog=None,FileName="Log",directory=".",ForceLog=False):
	if LogOn or ForceLog:
		if texttolog:
			buf=texttolog+"%s"%os.linesep
		else:
			buf=DescribeSituation()+"%s"%os.linesep
		f=open("%s%s%s.txt"%(directory,os.sep,FileName),'a')
		f.write(buf)
		f.close()
	


def combinations(iterable, r):
	# from the itertools for 2.6.1
	# combinations('ABCD', 2) --> AB AC AD BC BD CD
	# combinations(range(4), 3) --> 012 013 023 123
	if len(iterable)<r : return#I added this line. Pietro Speroni
	pool = tuple(iterable)
	n = len(pool)
	indices = range(r)
	yield tuple(pool[i] for i in indices)
	while 1:
		for i in reversed(range(r)):
			if indices[i] != i + n - r:
				break
		else:
			return
		indices[i] += 1
		for j in range(i+1, r):
			indices[j] = indices[j-1] + 1
		yield tuple(pool[i] for i in indices)


def findsubsets(S,m):#given the set S find all the subsets of order m
	return set(combinations(S, m))


###########################################
######## Database #########################
###########################################

def initalisedatabases(DirectoryAC):
	"""Call this when you open the databases for the first time"""
    #	global moleculelibrary, OrgLibrary, urllib, reactionDB, mutationsDB, LatticeLibrary, GeneratorsLibrary, OrganizationName, NameOrganization, TimeBiggerOrganizationName, IntersectionRelOrg, DirectoryAC
	global moleculelibrary, OrgLibrary, reactionDB, mutationsDB, LatticeLibrary, GeneratorsLibrary, OrganizationName, NameOrganization, TimeBiggerOrganizationName, IntersectionRelOrg
    #	DirectoryAC="%s%s%s%s%s"%(directory,os.sep,os.pardir,os.sep,urllib)
    #	if not os.path.isdir(DirectoryAC):
    #		os.mkdir(DirectoryAC)
    #if not os.path.isdir("%s%s%s%s%s"%(directory,os.sep,os.pardir,os.sep,urllib)):
	#		os.mkdir("%s%s%s%s%s"%(directory,os.sep,os.pardir,os.sep,urllib))
	moleculelibrary=MoleculeLibraryDB(DirectoryAC)
    #	moleculelibrary=MoleculeLibraryDB("%s%s%s%s%s"%(directory,os.sep,os.pardir,os.sep,urllib))
	OrgLibrary=OrganizationLibraryDB(DirectoryAC)
	reactionDB=MoleculeReactionsDB(DirectoryAC)   #this is not in the library as it depends on the AC. Note the Grammar error!
	mutationsDB=MoleculeMutationsDB(DirectoryAC)  #this is not in the library as it depends on the AC. Also it is not used if you just look for the Lattice
	LatticeLibrary=LatticeLibraryDB(DirectoryAC)  #a library from the set to the lattices. Probably too big?
	GeneratorsLibrary=GeneratorsLibraryDB(DirectoryAC)  #from an organization to the set of generators. Not used
	OrganizationName=OrganizationNameDB(DirectoryAC)    #very important
	NameOrganization=NameOrganizationDB(DirectoryAC)    #very important
	TimeBiggerOrganizationName=TimeBiggerOrganizationNameDB(DirectoryAC)   #a library I don't even know what it does and can be probably deleted!
	IntersectionRelOrg=RelationOrganizationDB("Intersection",DirectoryAC)    #probably not used as not memory efficient? In any case takes an org and responds with a couple of orgs that produce it via the intersection
	print "I am opening the archives"




def closedatabases():
	global moleculelibrary, OrgLibrary, reactionDB, mutationsDB, LatticeLibrary,TimeBiggerOrganizationName,IntersectionRelOrg
	moleculelibrary.DB.close()
	OrgLibrary.DB.close()
	reactionDB.DB.close()
	mutationsDB.DB.close()
	LatticeLibrary.DB.close()
	GeneratorsLibrary.DB.close()
	OrganizationName.DB.close()
	NameOrganization.DB.close()
	TimeBiggerOrganizationName.DB.close()
		
	IntersectionRelOrg.DB.close()
		
	global AllStates, AllRndSeeds
	AllStates.close()
	AllRndSeeds.close()
	print "I am closing the archives"
	

class MoleculeLibraryDB(object):
	"""a class the handles the Database of molecular reactions. Every time a reaction is tried it will call the function check.
	if the reaction is present in the DB, the result is returned. If not it is calculated, added to the DB and returned
	"""
	def __repr__(self):	 return 'MoleculeLibraryDB('+`self.DB`+')'
	def __str__(self):	 return self.DB.__str__()
	#	 def __init__(self,DB={}): self.DB=DB
	def __init__(self,directory=".",DB={}):
		self.DB=shelve.open("%s%sMoleculeLibrary.shelve"%(directory,os.sep),writeback=True)
		for k in DB.keys():
			self.DB[str(k)]=DB[k]
	def check(self,m0):
		try:
			return self.DB[str(m0)]
		except KeyError:
			molec=Molecule(m0)
			self.DB[str(m0)]=molec	#slightly faster by taking away a hash call
			return molec

class OrganizationLibraryDB(object):
	"""a class the handles the Database from set of molecules to organizations"""
	def __repr__(self):	 return 'OrganizationLibraryDB('+`self.DB`+')'
	def __str__(self):	 return self.DB.__str__()
	#	 def __init__(self,DB={}): self.DB=DB
	def __init__(self,directory=".",DB={}):
		self.DB=shelve.open("%s%sOrganizationLibrary.shelve"%(directory,os.sep),writeback=True)
		for k in DB.keys():
			self.DB[str(k)]=DB[k]
	def check(self,molecules):#I return the organization
		mollist=list(molecules)
		mollist.sort()
		frozm=str(mollist)
		try:			
			return NameOrganization.check(self.DB[frozm])	 #I return the organization
		except KeyError:
			Org=frozenset(Organization(molecules))
			name=OrganizationName.check(Org)
			self.DB[frozm]=name	 #slightly faster by taking away a hash call
			GeneratorsLibrary.add(name,molecules)
			return Org#I return the organization
	def checkname(self,molecules):#I return the name of the organization
		mollist=list(molecules)
		mollist.sort()
		frozm=str(mollist)
		#		 frozm=str(molecules)
		try:
			return self.DB[frozm]#I return the name of the organization
		except KeyError:
			Org=frozenset(Organization(molecules))
			name=OrganizationNameDB.check(Org)
			self.DB[frozm]=name	 #slightly faster by taking away a hash call
			GeneratorsLibrary.add(name,molecules)
			return name#I return the name of the organization


#TODO is this really necessary or is it too big and never called?
class LatticeLibraryDB(object):
	"""a class the handles the Database from set of molecules to organizations"""
	def __repr__(self):	 return 'OrganizationLibraryDB('+`self.DB`+')'
	def __str__(self):	 return self.DB.__str__()
	def __init__(self,DB={}): self.DB=DB
	def __init__(self,directory=".",DB={}):
		self.DB=shelve.open("%s%sLatticeLibrary.shelve"%(directory,os.sep),writeback=True)
		for k in DB.keys():
			self.DB[str(k)]=DB[k]
	def check(self,molecules):
		molecules=OrgLibrary.check(molecules)#this only tests on the biggest organisation
		try:
			result=self.DB[str(molecules)]
			return result
		except KeyError:
		#	 try:
		 #		 Os=CheckAllSets(frozenset(molecules))
		  #		 self.DB[str(molecules)]=Os
		   #	 print "added",len(Os)," orgs under",OrgLibrary.check(molecules)
			try:
				Os=CheckAllSets(frozenset(molecules))
				self.DB[str(molecules)]=Os
				print "added",len(Os)," orgs under",OrgLibrary.check(molecules)
			except OverflowError:
				Os=FindOrgKnownBelow(frozenset(OrgLibrary.check(molecules)))
			return Os

#TODO Also this is probably not used
class GeneratorsLibraryDB(object):
	"""a class the handles the Database from organization to the set of sets of molecules that generate them"""
	def __repr__(self):	 return 'GeneratorsLibraryDB('+`self.DB`+')'
	def __str__(self):	 return self.DB.__str__()
	def __init__(self,directory=".",DB={}):
		self.DB=shelve.open("%s%sGeneratorsLibrary.shelve"%(directory,os.sep),writeback=True)
		for k in DB.keys():
			self.DB[str(k)]=DB[k]
	def check(self,name):
		try:
			return self.DB[name]
		except KeyError:
			molecules=NameOrganization.check(name)
			generators=frozenset([molecules])
			self.DB[name]=generators  #slightly faster by taking away a hash call
			return generators
	def add(self,name,molecules):
		name=str(name)
		molecules=frozenset(molecules)
		try:
			generators=self.DB[name]|frozenset([molecules])
			self.DB[name]=generators 
			return generators
		except KeyError:
			O=NameOrganization.check(name)
			generators=frozenset([O])|frozenset([molecules])
			self.DB[name]=generators 
			return generators

class OrganizationNameDB(object):
	"""a class that handles the relation Name-Organization, takes a frozenset and returns a string"""
	#	 def __repr__(self):  return 'OrganizationLibraryDB('+`self.DB`+')'
	def __repr__(self):	 return 'OrganizationNameDB('+`self.DB`+')'
	def __str__(self):	 return self.DB.__str__()
	def __init__(self,directory=".",DB={}):
		self.DB=shelve.open("%s%sOrganizationName.shelve"%(directory,os.sep),writeback=True)
		for k in DB.keys():
			self.DB[str(k)]=DB[k]
			#		 OrganizationNameDB.max="aaaaa"
		OrganizationNameDB.max=0
		values=self.DB.values()
		if len(values):
			values=[int(v) for v in values]
			OrganizationNameDB.max=max(values)
		print OrganizationNameDB.max,
	def backtrack(self,name):
		org=[k for k, v in self.DB.iteritems() if v == name]
		TempOrg="org="+org[0]
		exec TempOrg
		org=frozenset(org)
		return org
		#return org[0]
	def check(self,org):
		orglist=list(org)
		orglist.sort()
		orgstr=str(orglist)
		try:
			return int(self.DB[orgstr])
		except KeyError:
			#			 OrganizationNameDB.max=next_string(OrganizationNameDB.max)
			OrganizationNameDB.max+=1
			self.DB[orgstr]=str(OrganizationNameDB.max)
			NameOrganization.add(OrganizationNameDB.max,org)
			return OrganizationNameDB.max


class NameOrganizationDB(object):
	"""a class the handles the relation Name-Organization, takes a string and returns a frozenset"""
	def __repr__(self):	 return 'NameOrganizationDB('+`self.DB`+')'
	#	 def __repr__(self):  return 'OrganizationLibraryDB('+`self.DB`+')'
	def __str__(self):	 return self.DB.__str__()
	def __init__(self,directory=".",DB={}):
		self.DB=shelve.open("%s%sNameOrganization.shelve"%(directory,os.sep),writeback=True)
		for k in DB.keys():
			self.DB[str(k)]=DB[k]
	def check(self,name):
		name=str(name)
		try:
			org=self.DB[name]
			if type(org)==type(""):
				print "I stored it as a string, let me change that"
				Temp="org="+org
				exec(Temp)
				org=frozenset(Temp)
				self.DB[name]=org
			elif type(org)==type(set([])):
				print "I stored it as a set, let me change that"
				org=frozenset(org)
				self.DB[name]=org
			return org
			#			return self.DB[name]
		except KeyError:
			org=OrganizationName.backtrack(name)
			self.DB[name]=org
			print "aaargh ERROR 783942,", name, org
			return org
	def add(self,name,org):
		if type(org)==type(set([])):
			org=frozenset(org)
		assert type(org)==type(frozenset([]))
		self.DB[str(name)]=org

#TODO Also this is probably not used
class TimeBiggerOrganizationNameDB(object):
	"""a class the handles the relation time to the Organization (by name), takes a number and returns the string of the organization"""
	def __repr__(self):	 return 'TimeBiggerOrganizationNameDB('+`self.DB`+')'
	def __str__(self):	 return self.DB.__str__()
	def __init__(self,directory=".",DB={}):
		self.DB=shelve.open("%s%sTimeBiggerOrganizationName.shelve"%(directory,os.sep),writeback=True)
		for k in DB.keys():
			self.DB[str(k)]=DB[k]
	def check(self,step):
		try:
			return self.DB[step]
		except KeyError:
			#			 print "aaargh ERROR 23243"
			b=Cell(AllStates["%10s"%step])
			self.DB[step]=OrganizationLibrary.checkname(b.molecules)
	def add(self,time,orgName):
		self.DB[step]=orgName




#TODO Also this is probably not used, but at the moment is only called for the intersection, meaning you give an organisations and will tell you couples of organisations whose intersection is that organisation
class RelationOrganizationDB(object):
	"""a class the handles the relation set of organizations, organization, takes a frozenset and returns an organization"""
	def __repr__(self):	 return 'RelationOrganizationDB('+`self.DB`+')'
	def __str__(self):	 return self.DB.__str__()
	def __init__(self,nameRelation,directory=".",DB={}):
		self.DB=shelve.open("%s%s%sRelOrg.shelve"%(directory,os.sep,nameRelation),writeback=True)
		for k in DB.keys():
			self.DB[str(k)]=DB[k]
	def check(self,name):
		name=str(name)
		try:
			return self.DB[name]
		except KeyError:
			return "none"
	def add(self,name,org):
		self.DB[str(name)]=org
	def backtrack(self,name):
		"""given an organizations tells you a series of couples of organizations that returns that"""
		#keysinteresting=find_keys(self.DB, name)
		keysinteresting=[k for k, v in self.DB.iteritems() if v == name and k[0]!= name and k[1]!=name]
		keyinter=keysinteresting[0]
		OrgA=NameOrganization.check(keyinter[0])
		OrgB=NameOrganization.check(keyinter[1])
		OrgAiOrgB = OrgLibrary.check(OrgA & OrgB)
		return OrgAiOrgB		





def setuparchives(directory):
	""" Call this before you start, to store all the states"""
	global AllStates, AllRndSeeds
	print "I am setting up the ALLSTATES archive"
	#	AllStates=shelve.open("%s%sAllStates.shelve"%(directory,os.sep),writeback=False)     
	#	AllStates=shelve.open("%s%sAllStates1.shelve"%(directory,os.sep))     
	AllStates=shelve.open("%s%sAllStates.shelve"%(directory,os.sep))     
	#if writeback=True the DB gets corrupted when I access all the entries
	AllRndSeeds=shelve.open("%s%sAllRndSeeds.shelve"%(directory,os.sep),writeback=True)


###########################################
######## Reaction Network #################
###########################################

def ReactSets(molLeft,molRight):
	"""Returns the list of all the molecules generated by a set of molecules, with their relative multiplicity
	"""
	results=set([])
	for m in molLeft:
		for n in molRight:
			res=react(m,n)
			if res!= None:
				results|=set([res])
	return results



###########################################
######## Chemical Organization Theory #######
###########################################

def Closure(molecules):
	"""
	Returns the set of the closure of a given list of molecules
	"""
	newmol=set(molecules)
	oldmol=set([])
	while newmol:
		gen=ReactSets(newmol,newmol)
		gen|=ReactSets(newmol,oldmol)
		gen|=ReactSets(oldmol,newmol)
		oldmol|=newmol
		newmol=gen-oldmol
	return oldmol

def SelfMaintainance(molecules):
	"""
	Returns the biggest self maintaining set inside a given list of molecules
	"""
	mol=set(molecules)
	l=len(mol)
	while l:
		notgeneratedmolecules=copy.deepcopy(mol)
		for t in mol:
			for s in mol:
				r=react(t,s)
				if r!= None:
					notgeneratedmolecules-=set([r])
					if not len(notgeneratedmolecules):
						return mol
		mol-=notgeneratedmolecules
		l=len(mol)
	return set([])

def Organization(molecules):
	"""
	Returns the organization generated by a set of molecules
	"""
	return SelfMaintainance(Closure(molecules))


def IsOrganization(molecules):
	"""
	Returns True if the set of molecules is an organization. It is about one order of magnitude faster (but messiers) than the previous one as it checks while it is running.
	"""
	mol=set(molecules)
	tobegenerated=copy.deepcopy(mol)
	for m in mol:
		for n in mol:
			res=react(m,n)
			if res!= None:
				if res not in mol:
					return False #Not Closed
				tobegenerated -= set([res])
	if tobegenerated:
		return False #Not Self Mainaining
	return True


###########################################
######## Calculate Organization #######
###########################################

def IsOrg1inOrg2Simple(A,B):
	"""
	Given two organizations, A and B, returns true if A < B,
	This is equivalent to check if (A intersection B) = A
	COST: SIZE A+SIZE B,+SET COMPARISON of SIZE A 
	
	"""
	OrgA=NameOrganization.check(A)
	OrgB=NameOrganization.check(B)
	if (OrgA & OrgB)==OrgA:
		return 1
	return 0


def Org1IntersectionOrg2(OrgA,OrgB):
	OrgDetailsA=NameOrganization.check(OrgA)
	OrgDetailsB=NameOrganization.check(OrgB)
	OrgResultDetails=Organization(OrgDetailsA&OrgDetailsB)
	OrgResult=OrganizationName.check(OrgResultDetails)
	return OrgResult

def Org1UnionOrg2(OrgA,OrgB):
	OrgDetailsA=NameOrganization.check(OrgA)
	OrgDetailsB=NameOrganization.check(OrgB)
	OrgResultDetails=Organization(OrgDetailsA|OrgDetailsB)
	OrgResult=OrganizationName.check(OrgResultDetails)
	return OrgResult








###########################################
######## Map Organizations Brute Force#######
###########################################


def powerset(iterable):
    """powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"""
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))

def subsets(s):
    return map(set, powerset(s))


def CheckAllSubsets(Molecules):
	"Checks one by one all the subsets of molecules, to find which one are organisations. But then returns them without storing them"
	organisationsFound=set([])
	subsetsMol=subsets(Molecules)
	for s in subsetsMol:
		if IsOrganization(s):
			organisationsFound.add(frozenset(s))
			print "Organisation found", len(organisationsFound), s
	return organisationsFound
		

###########################################
######## Map Organizations #######
###########################################

NameSize={}  #often it is very easy if we can simply associate to an org the size, this shortcuts many calculations


#We need to store, for every organization, the organizations that are above, below, strictly above (covered), and so on.
OCovered={}
OMaybeCovered={}
OCovering={}
OMaybeCovering={}
OBelowNotCovered={}
OSide={}
OAboveNotCovering={}

OrganizationsStudied=set([])

TableIntersection=SymetricTable() #Given two organizations (A,B), stored by name, returns the union of them
TableUnion=SymetricTable() #Given two organizations (A,B), stored by name, returns the union of them. Stores each instance just once unsig the fact that is symmetric


#given an organization, we need to map if we add molecules, where would we be going
DownwardMolecules=dict()
UpwardMolecules=dict()
SidewardMolecules=dict()



UnionTriangulated=0
UnionByAUBBuilding=0
UnionByAUBUC1=0
UnionByAUBUC2=0
UnionByAUBcovered=0
UnionFoundByHand=0
NewOrgUnionFoundByHand=0

IntersectionByAIBBuilding=0
IntersectionTriangulated=0
IntersectionByAIntBIntC1=0
IntersectionByAIntBIntC2=0
IntersectionByAIBcovered=0
IntersectionFoundByHand=0
NewOrgIntersectionFoundByHand=0

ToBeSolvedByHandUnionsDict={}
ToBeSolvedByHandIntersectionsDict={}


def DescribeSituation():
	"""Describe Situation returns a string describing the number of organizations, how many in the tables union and intersection, and how many were calculated by hand versus how many were found through theorems.
	"""
	buf=""
	
	
	len_OrganizationsStudied=len(OrganizationsStudied)
	SizeTables=len_OrganizationsStudied+len_OrganizationsStudied*(len_OrganizationsStudied-1)/2
	#	buf+="TI= %s  "%len(TableIntersection)
	buf+="TI=%.0f%% "%((100*len(TableIntersection))/SizeTables)
	#	buf+="TU= %s  "%len(TableUnion)
	buf+="TU=%.0f%% "%((100*len(TableUnion))/SizeTables)
	#	buf+="STU= %s  "%len(STableUnion)
	
	buf+="#O=%s "%len(OrganizationsStudied)
	AlgebraUnionResults=float(UnionTriangulated+UnionByAUBUC1+UnionByAUBUC2+UnionByAUBcovered+UnionByAUBBuilding)
	HandUnionResults=float(UnionFoundByHand)
	NewHandUnionResults=float(NewOrgUnionFoundByHand)
	AllHandUnionResults=float(HandUnionResults+NewHandUnionResults)
	
	TotalUnionResults=float(AlgebraUnionResults+AllHandUnionResults)
	#	TotalUnionResults=AlgebraUnionResults+HandUnionResults+NewHandUnionResults
	#	if TotalUnionResults:		buf+=" U(A:H:N)= %.2f%% %.2f%% %.2f%% "%((AlgebraUnionResults/TotalUnionResults)*100,(HandUnionResults/TotalUnionResults)*100,(NewHandUnionResults/TotalUnionResults)*100)
	if TotalUnionResults:
		#		buf+=" U(AB:AT:AU1:AU2:AC:H:N)= %.2f%% %.2f%% %.2f%% %.2f%% %.2f%% %.2f%% %.2f%% "%((UnionByAUBBuilding/TotalUnionResults)*100,(UnionTriangulated/TotalUnionResults)*100,(UnionByAUBUC1/TotalUnionResults)*100,(UnionByAUBUC2/TotalUnionResults)*100,(UnionByAUBcovered/TotalUnionResults)*100,(HandUnionResults/TotalUnionResults)*100,(NewHandUnionResults/TotalUnionResults)*100)
		buf+=" U(A:H)= %.2f%% %.2f%% "%((AlgebraUnionResults/TotalUnionResults)*100,(AllHandUnionResults/TotalUnionResults)*100)
		
	AlgebraIntersectionResults=float(IntersectionTriangulated+IntersectionByAIntBIntC1+IntersectionByAIntBIntC2+IntersectionByAIBcovered+IntersectionByAIBBuilding)
	HandIntersectionResults=float(IntersectionFoundByHand)
	NewHandIntersectionResults=float(NewOrgIntersectionFoundByHand)
	AllHandUnionResults=float(HandIntersectionResults+NewHandIntersectionResults)
	TotalIntersectionResults=AlgebraIntersectionResults+HandIntersectionResults+NewHandIntersectionResults
	if TotalIntersectionResults:
		#buf+="^(A:H:N)= %.2f%% %.2f%% %.2f%% "%((AlgebraIntersectionResults/TotalIntersectionResults)*100,(HandIntersectionResults/TotalIntersectionResults)*100,(NewHandIntersectionResults/TotalIntersectionResults)*100)
		#buf+=" ^(AB:AT:AI1:AI2:AC:H:N)= %.2f%% %.2f%% %.2f%% %.2f%% %.2f%% %.2f%% %.2f%% "%((IntersectionByAIBBuilding/TotalIntersectionResults)*100,(IntersectionTriangulated/TotalIntersectionResults)*100,(IntersectionByAIntBIntC1/TotalIntersectionResults)*100,(IntersectionByAIntBIntC2/TotalIntersectionResults)*100,(IntersectionByAIBcovered/TotalIntersectionResults)*100,(HandIntersectionResults/TotalIntersectionResults)*100,(NewHandIntersectionResults/TotalIntersectionResults)*100)
		buf+=" ^(A:H)= %.2f%% %.2f%% "%((AlgebraIntersectionResults/TotalIntersectionResults)*100,(AllHandUnionResults/TotalIntersectionResults)*100)
	
	global starttime
	averagetime=(time.clock()-starttime)/(len_OrganizationsStudied-1)
	buf+=" time= %s averagetime=%s "%(time.clock()-starttime, averagetime)
		
	#	buf+=" time= %s d %s h %s' %s\""%(time.localtime().tm_yday, time.localtime().tm_hour, time.localtime().tm_min, time.localtime().tm_sec)
	return buf
	
	
	
def SolveByHandUnion(Org):
	"""sometimes you just need to calculate the union by hand"""
	global ToBeSolvedByHandIntersectionsDict
	global ToBeSolvedByHandUnionDict
	NewOrganizations=set([])
					
	#buf=TestAllRelations();	
	#if buf:			print buf
	#else: print "all relations are fine."
	
	Orgs=ToBeSolvedByHandUnionsDict.pop(Org)
	OrgsSizePairs=[(o,NameSize[o]) for o in Orgs]
	OrgsSizePairs.sort(sortpairs, reverse=False)	
	for OrgP,s in OrgsSizePairs:
		if (Org,OrgP) in TableUnion:	continue
		if AddOneOrgUnion(Org,OrgP):	
			#sys.stdout.write('u'); sys.stdout.flush()		
			continue
		if AddOneOrgUnion(OrgP,Org):	
			#sys.stdout.write('v');sys.stdout.flush()		
			continue
		OrgC=CheckUnion(Org,OrgP)
		if OrgC:
			print DescribeSituation(),
			sys.stdout.flush()		
			ToBeSolvedByHandUnionsTemp,ToBeSolvedByHandIntersectionsTemp=AddOrganizationOfficiallyFromUnion(OrgC,Org,OrgP)
			if ToBeSolvedByHandUnionsTemp:					ToBeSolvedByHandUnionsDict[OrgC]=ToBeSolvedByHandUnionsTemp
			if ToBeSolvedByHandIntersectionsTemp:			ToBeSolvedByHandIntersectionsDict[OrgC]=ToBeSolvedByHandIntersectionsTemp
			NewOrganizations|=set([OrgC])
			AddUnionToTable(OrgP,Org,OrgC)
		else:
			#sys.stdout.write('U');sys.stdout.flush()		
			AddUnionToTable(OrgP,Org,TableUnion[(OrgP,Org)])
		sys.stdout.flush()
	return NewOrganizations

def SolveByHandIntersection(Org):
	"""sometimes you just need to calculate the intersection by hand"""
	
	global ToBeSolvedByHandIntersectionsDict
	global ToBeSolvedByHandUnionsDict
	NewOrganizations=set([])
	Orgs=ToBeSolvedByHandIntersectionsDict.pop(Org)
	OrgsSizePairs=[(o,NameSize[o]) for o in Orgs]
	OrgsSizePairs.sort(sortpairs, reverse=True)	
	for OrgP,s in OrgsSizePairs:
		if (Org,OrgP) in TableIntersection:		continue
		if AddOneOrgIntersection(Org,OrgP):		
			#sys.stdout.write('i') ;			sys.stdout.flush()		
			continue
		if AddOneOrgIntersection(OrgP,Org):		
			#sys.stdout.write('j') ;sys.stdout.flush()		
			continue			
		OrgC=CheckIntersection(Org,OrgP)
		if OrgC:
			print DescribeSituation(),
			sys.stdout.flush()
			ToBeSolvedByHandUnionsTemp,ToBeSolvedByHandIntersectionsTemp=AddOrganizationOfficiallyFromIntersection(OrgC,Org,OrgP)
			if ToBeSolvedByHandUnionsTemp:					ToBeSolvedByHandUnionsDict[OrgC]=ToBeSolvedByHandUnionsTemp
			if ToBeSolvedByHandIntersectionsTemp:			ToBeSolvedByHandIntersectionsDict[OrgC]=ToBeSolvedByHandIntersectionsTemp
			NewOrganizations|=set([OrgC])	
			AddIntersectionToTable(OrgP,Org,OrgC)
		else:
			#sys.stdout.write('I');			sys.stdout.flush()		
			AddIntersectionToTable(OrgP,Org,TableIntersection[(OrgP,Org)])		
	sys.stdout.flush()			
	return NewOrganizations
	
	
	
	
def CheckUnion(OrgA,OrgB):
	"""Looks for the union of two organizations by calculating it and then adding it"""
	if (OrgA, OrgB) in TableUnion:	return None
	global 	OrganizationsStudied
	global	NewOrgUnionFoundByHand
	global	UnionFoundByHand
	OrgC=Org1UnionOrg2(OrgA,OrgB)
	if OrgC in OrganizationsStudied:#		print "Organization ", OrgC, "found by hand as ",OrgA,"union",OrgB
		UnionFoundByHand+=1
		TableUnion[(OrgA,OrgB)]=OrgC
		#		TableUnion[(OrgB,OrgA)]=OrgC
		return None
	else:
		print "New Org: ", OrgC, " := ",OrgA," U ",OrgB
		TableUnion[(OrgA,OrgB)]=OrgC
		#		TableUnion[(OrgB,OrgA)]=OrgC	
		NewOrgUnionFoundByHand+=1		
		return OrgC


def CheckIntersection(OrgA,OrgB):
	"""Looks for the Intersection of two organizations by calculating it and then adding it"""
	global 	OrganizationsStudied
	global	NewOrgIntersectionFoundByHand
	global	IntersectionFoundByHand
	if (OrgA, OrgB) in TableIntersection: #		print "I dont' need to look for the new org"
		return None
	OrgC=Org1IntersectionOrg2(OrgA,OrgB)
	if OrgC in OrganizationsStudied:#		print "Organization ", OrgC, "found by hand as ",OrgA,"intersection",OrgB
		IntersectionFoundByHand+=1
		TableIntersection[(OrgA,OrgB)]=OrgC
		#		TableIntersection[(OrgB,OrgA)]=OrgC
		return None
	else:
		print "New Org: ", OrgC, " := ",OrgA," ^ ",OrgB
		TableIntersection[(OrgA,OrgB)]=OrgC
		#		TableIntersection[(OrgB,OrgA)]=OrgC	
		NewOrgIntersectionFoundByHand+=1		
		return OrgC





def FindGenerators(Orgs,Relation,ToBeGenerated):
	"""This powerful function takes a set of organizations, a relation (in the form of a dictionary), and another organization, and returns all the pair of organization from the 
	set such that the relation applied to the pair returns the wanted organization The general form to call it is:
	if len(Orgs)>1:
		for t in FindGenerators(Orgs,TableUnion,Org):
			print t
	"""
	#	if len(Orgs)<=1: return
	for s in (t for t in itertools.combinations(Orgs, 2) if Relation.get(t, None)==ToBeGenerated):
		yield s

def DiscoverGenerators(OrgsNear,OrgsFar):
	"""This function takes a couple of sets, and returns the pair with the first element coming from the first set, the second from the union of the first and the second
	but such that the two elements are never the same. It is used to look for pair of organizations that generate another organization
	"""
	AllOrgs=OrgsFar|OrgsNear
	for n in OrgsNear:
		AllOrgs-=set([n])
		for m in AllOrgs:
			yield (n,m)

	
	
def AddOneOrgIntersection(Org,OrgP): 
	"""Using only Algebraic means can you find the Intersection of Org and OrgP?
	takes Org and OrgP, and returns the result, if found, and the unsolved intersection cases it might have met in trying to integrate those elements
	"""
	try:					return TableIntersection[(Org,OrgP)]
	except KeyError:		pass
	global IntersectionTriangulated
	global IntersectionByAIntBIntC1
	global IntersectionByAIntBIntC2
	global IntersectionFoundByHand
	global NewOrgIntersectionFoundByHand
	#	OrgsIntersectionUpperApprox=set([])
	#	for OrgQ in OCovering[Org]|OMaybeCovering[Org]:
	#		try:				OrgsIntersectionUpperApprox|=set([TableIntersection[(OrgP,OrgQ)]])
	#		except KeyError:	pass
	#	if OrgsIntersectionUpperApprox:
	#		for OrgQ in OCovered[Org]|OMaybeCovered[Org]:
	#			try:				OrgIntersectionLowerApprox=TableIntersection[(OrgP,OrgQ)]
	#			except KeyError:	continue
	#			if OrgIntersectionLowerApprox in OrgsIntersectionUpperApprox:
	#				TableIntersection[(OrgP,Org)]=OrgIntersectionLowerApprox
	#				TableIntersection[(Org,OrgP)]=OrgIntersectionLowerApprox
	#				IntersectionTriangulated+=1
	#				return OrgIntersectionLowerApprox

	OrgsAboveOrg			 = OrganizationsAbove(Org)
	OrgsAboveOrgP			 = OrganizationsAbove(OrgP)
	#	OrgsAboveOrg             = OCovering[Org]  | OMaybeCovering[Org]
	#	OrgsAboveOrgP            = OCovering[OrgP] | OMaybeCovering[OrgP]
	OrgsAboveOrgInteresting  = OrgsAboveOrg.difference(OrgsAboveOrgP)
	OrgsAboveOrgPInteresting = OrgsAboveOrgP.difference(OrgsAboveOrg)

	OrgsBelowOrg             = OrganizationsUnder(Org)
	OrgsBelowOrgP            = OrganizationsUnder(OrgP)
	#	OrgsBelowOrgP            = OCovered[OrgP]  | OMaybeCovered[OrgP]
	#	OrgsBelowOrg             = OCovered[Org]   | OMaybeCovered[Org]
	OrgsBelowOrgP            = OCovered[OrgP]  | OMaybeCovered[OrgP]
	OrgsBelowOrgInteresting  = OrgsBelowOrg.difference(OrgsBelowOrgP)
	OrgsBelowOrgPInteresting = OrgsBelowOrgP.difference(OrgsBelowOrg)

	OrgsToIgnore = set([])
	for OrgQ in OrgsAboveOrgInteresting:    #I take all the organizations Above Org
		if OrgQ in OrgsToIgnore:
			continue
		try:
			OrgIntersection=TableIntersection[(OrgP,OrgQ)] 			
			if OrgIntersection in OrgsBelowOrg:
				TableIntersection[(OrgP,Org)]=OrgIntersection
			#	TableIntersection[(Org,OrgP)]=OrgIntersection
				if ShowIntersection:
					#					print "%s ^ %s, but %s > %s and %s ^ %s = %s, with %s < %s thus %s ^ %s = %s "%(Org, OrgP, OrgQ, Org, OrgQ, OrgP, OrgIntersection, OrgIntersection, Org, Org, OrgP, OrgIntersection)
				    print "%s ^ %s, but %s > %s > %s and %s ^ %s = %s, thus %s ^ %s = %s "%(Org, OrgP, OrgQ, Org, OrgIntersection, OrgQ, OrgP, OrgIntersection, Org, OrgP, OrgIntersection)					
				IntersectionTriangulated+=1
				return OrgIntersection
		except KeyError:	
			OrgsToIgnore |= OrganizationsAbove(OrgQ)
			continue				
			
	OrgsToIgnore = set([])
	for OrgQ in OrgsAboveOrgPInteresting:    #I take all the organizations Above Org
		if OrgQ in OrgsToIgnore:
			continue
		try:
			OrgIntersection=TableIntersection[(Org,OrgQ)]
			if OrgIntersection in 	OrgsBelowOrgP:
				TableIntersection[(OrgP,Org)]=OrgIntersection
			#	TableIntersection[(Org,OrgP)]=OrgIntersection
				if ShowIntersection:
				    print "%s ^ %s, but %s > %s > %s and %s ^ %s = %s, thus %s ^ %s = %s "%(OrgP, Org, OrgQ, OrgP, OrgIntersection, OrgQ, Org, OrgIntersection, OrgP, Org, OrgIntersection)					
				IntersectionTriangulated+=1
				return OrgIntersection
		except KeyError:	
			OrgsToIgnore |= OrganizationsAbove(OrgQ)
			continue				

	for (OrgR,OrgS) in DiscoverGenerators(OCovering[Org],OMaybeCovering[Org]):
			TableIntersection[(OrgR,OrgS)]=Org
	#		TableIntersection[(OrgS,OrgR)]=Org			
			try:
				try:					OrgS_inter_OrgP=TableIntersection[(OrgS,OrgP)]									
				except KeyError:		raise
				try:					ResOrg=TableIntersection[(OrgR,OrgS_inter_OrgP)]
				except KeyError:		raise
				TableIntersection[(OrgP,Org)]=ResOrg
			#	TableIntersection[(Org,OrgP)]=ResOrg
				IntersectionByAIntBIntC1+=1
				if ShowIntersection:
					print "IntersectionByAIntBIntC1"				
				return ResOrg
			except KeyError:													
				#		(OrgS ^ OrgR) ^ OrgP    =   OrgS ^ (OrgR ^ OrgP);but Organization=OrgR ^ OrgS;thus Organization ^ OrgP=OrgR ^ (OrgS ^ OrgP)
				try:					OrgR_inter_OrgP=TableIntersection[(OrgR,OrgP)]									
				except KeyError:		continue
				try:					ResOrg=TableIntersection[(OrgS,OrgR_inter_OrgP)]
				except KeyError:		continue
				TableIntersection[(OrgP,Org)]=ResOrg
			#	TableIntersection[(Org,OrgP)]=ResOrg
				IntersectionByAIntBIntC2+=1
				if ShowIntersection:
					print "IntersectionByAIntBIntC2"
				return ResOrg
	if len(OMaybeCovering[Org])>1:
		for (OrgS,OrgR) in FindGenerators(OMaybeCovering[Org],TableIntersection,Org):
			try:
				try:				OrgS_inter_OrgP=TableIntersection[(OrgS,OrgP)]									
				except KeyError:	raise
				try:				ResOrg=TableIntersection[(OrgR,OrgS_inter_OrgP)]
				except KeyError:	raise
				TableIntersection[(OrgP,Org)]=ResOrg
			#	TableIntersection[(Org,OrgP)]=ResOrg
				IntersectionByAIntBIntC1+=1
				if ShowIntersection:
					print "IntersectionByAIntBIntC1"
				return ResOrg
			except KeyError:													#		(OrgS ^ OrgR) ^ OrgP    =   OrgS ^ (OrgR ^ OrgP);but Organization=OrgR ^ OrgS;thus Organization ^ OrgP=OrgR ^ (OrgS ^ OrgP)
				try:				OrgR_inter_OrgP=TableIntersection[(OrgR,OrgP)]									
				except KeyError:	continue
				try:				ResOrg=TableIntersection[(OrgS,OrgR_inter_OrgP)]
				except KeyError:	continue
				TableIntersection[(OrgP,Org)]=ResOrg
		#		TableIntersection[(Org,OrgP)]=ResOrg
				IntersectionByAIntBIntC2+=1
				if ShowIntersection:
					print "IntersectionByAIntBIntC2"
				return ResOrg
	# If the lattice was distributive here you could add that if A, B are such that A U B = Org then (A ^ OrgP) U (B ^ OrgP)=(A U B)^ OrgP= Org ^ OrgP
	return None

def AddOneOrgUnion(Org,OrgP):
	"""Using only Algebraic means can you find the Union of Org and OrgP?
	takes Org and OrgP, and returns the result, if found, and the unsolved Union cases it might have met in trying to integrate those elements
	"""
	try:				return TableUnion[(Org,OrgP)]	
	except KeyError:	pass
	#	try:				return STableUnion[(Org,OrgP)]	
	#	except KeyError:	pass
	global	UnionTriangulated
	global	UnionByAUBUC1
	global	UnionByAUBUC2
	global	UnionFoundByHand
	global  UnionByAUBcovered
	OrgsUnionLowApprox=set([])
	#	for OrgQ in OCovered[Org]|OMaybeCovered[Org]:    #I take all the organizations BELOW Org
	#		try:					OrgsUnionLowApprox|=set([TableUnion[(OrgP,OrgQ)]])   #For each of those organizations OrgQ I check if the union with the OrgP is above Org
	#		except KeyError:		pass
	#	if OrgsUnionLowApprox:
	#		for OrgQ in OCovering[Org]|OMaybeCovering[Org]:
	#			try:				OrgUnionUpperApprox=TableUnion[(OrgP,OrgQ)]
	#			except KeyError:	continue
	#			if OrgUnionUpperApprox in OrgsUnionLowApprox:
	#				TableUnion[(OrgP,Org)]=OrgUnionUpperApprox
	#				TableUnion[(Org,OrgP)]=OrgUnionUpperApprox
	#				if ShowUnion:
	#				    print "UnionTriangulated: %s U %s, but exits OrgQ such that  = %s"%(Org,OrgP,OrgUnionUpperApprox)
	#				UnionTriangulated+=1							
	#				return OrgUnionUpperApprox
				
				
				#if I find an organization below org that with the union with orgP makes something above org, then the union is the same
	OrgsAboveOrg			 = OrganizationsAbove(Org)
	OrgsAboveOrgP			 = OrganizationsAbove(OrgP)
	#	OrgsAboveOrg             = OCovering[Org]  | OMaybeCovering[Org]
	#	OrgsAboveOrgP            = OCovering[OrgP] | OMaybeCovering[OrgP]
	OrgsAboveOrgInteresting  = OrgsAboveOrg.difference(OrgsAboveOrgP)
	OrgsAboveOrgPInteresting = OrgsAboveOrgP.difference(OrgsAboveOrg)

	OrgsBelowOrg             = OrganizationsUnder(Org)
	OrgsBelowOrgP            = OrganizationsUnder(OrgP)
	#	OrgsBelowOrgP            = OCovered[OrgP]  | OMaybeCovered[OrgP]
	#	OrgsBelowOrg             = OCovered[Org]   | OMaybeCovered[Org]
	OrgsBelowOrgP            = OCovered[OrgP]  | OMaybeCovered[OrgP]
	OrgsBelowOrgInteresting  = OrgsBelowOrg.difference(OrgsBelowOrgP)
	OrgsBelowOrgPInteresting = OrgsBelowOrgP.difference(OrgsBelowOrg)

	OrgsToIgnore = set([])
	for OrgQ in OrgsBelowOrgInteresting:    #I take all the organizations BELOW Org
		if OrgQ in OrgsToIgnore:
			continue
		try:
			OrgUnion=TableUnion[(OrgP,OrgQ)]
			if OrgUnion in 	OrgsAboveOrg:
				TableUnion[(OrgP,Org)]=OrgUnion
				#				TableUnion[(Org,OrgP)]=OrgUnion
				#				STableUnion[(Org,OrgP)]=OrgUnion
				
				if ShowUnion:
				    print "%s U %s, but %s < %s < %s and %s U %s = %s, thus %s U %s = %s "%(OrgP, Org, OrgQ, Org, OrgUnion, OrgQ, OrgP, OrgUnion, OrgP, Org, OrgUnion)				
				UnionTriangulated+=1
				return OrgUnion
		except KeyError:	
			OrgsToIgnore |= OrganizationsUnder(OrgQ)
			continue				

	OrgsToIgnore = set([])
	for OrgQ in OrgsBelowOrgPInteresting:    #I take all the organizations BELOW Org
		if OrgQ in OrgsToIgnore:
			continue
		try:
			OrgUnion=TableUnion[(Org,OrgQ)]
			if OrgUnion in 	OrgsAboveOrgP:
				TableUnion[(OrgP,Org)]=OrgUnion
				#				TableUnion[(Org,OrgP)]=OrgUnion
				if ShowUnion:
				    print "%s U %s, but %s < %s < %s and %s U %s = %s, thus %s U %s = %s "%(Org, OrgP, OrgQ, OrgP, OrgUnion, OrgQ, Org, OrgUnion, Org, OrgP, OrgUnion)					
				UnionTriangulated+=1
				return OrgUnion
		except KeyError:	
			OrgsToIgnore |= OrganizationsUnder(OrgQ)
			continue				

				
	for (OrgR,OrgS) in DiscoverGenerators(OCovered[Org],OMaybeCovered[Org]):	
			try:				
				OrgRUOrgS=TableUnion[(OrgR,OrgS)]  #the union between two organizations, one that is contained in an organization and another one  
				assert (Org==OrgRUOrgS)									
			except KeyError:
				TableUnion[(OrgR,OrgS)]=Org   #because we take organization that are strictly contained in Org, the union should be = Org
				#				TableUnion[(OrgS,OrgR)]=Org
				UnionByAUBcovered+=1
			try:
				try:				OrgSUOrgP=TableUnion[(OrgS,OrgP)]									
				except KeyError:	raise
				try:				ResOrg=TableUnion[(OrgR,OrgSUOrgP)]
				except KeyError:	raise
				TableUnion[(OrgP,Org)]=ResOrg
				#				TableUnion[(Org,OrgP)]=ResOrg
				if ShowUnion:
				    print "%s U %s = (%s U %s) U %s = %s U (%s U %s) = %s U %s = %s"%(Org,OrgP,OrgR,OrgS,OrgP,OrgR,OrgS,OrgP,OrgR,OrgSUOrgP,ResOrg)
				UnionByAUBUC1+=1
				return ResOrg
			except KeyError:					
				try:	OrgRUOrgP=TableUnion[(OrgR,OrgP)]									
				except KeyError:										continue
				try:	ResOrg=TableUnion[(OrgS,OrgRUOrgP)]
				except KeyError:										continue			
				TableUnion[(OrgP,Org)]=ResOrg
				#				TableUnion[(Org,OrgP)]=ResOrg
				UnionByAUBUC2+=1
				if (ShowUnion): 
				    print "%s U %s = (%s U %s) U %s = %s U (%s U %s) = %s U %s = %s"%(Org,OrgP,OrgS,OrgR,OrgP,OrgS,OrgR,OrgP,OrgS,OrgRUOrgP,ResOrg)
				return ResOrg
	if len(OMaybeCovered[Org])>1:
		for (OrgS,OrgR) in FindGenerators(OMaybeCovered[Org],TableUnion,Org):
			try:
				try:				OrgSUOrgP=TableUnion[(OrgS,OrgP)]									
				except KeyError:	raise
				try:				ResOrg=TableUnion[(OrgR,OrgSUOrgP)]
				except KeyError:	raise
				TableUnion[(OrgP,Org)]=ResOrg
				#				TableUnion[(Org,OrgP)]=ResOrg
				if ShowUnion: print "%s U %s = (%s U %s) U %s = %s U (%s U %s) = %s U %s = %s"%(Org,OrgP,OrgR,OrgS,OrgP,OrgR,OrgS,OrgP,OrgR,OrgSUOrgP,ResOrg)
				UnionByAUBUC1+=1
				return ResOrg
			except KeyError:													#		(OrgS U OrgR) U OrgP    =   OrgS U (OrgR U OrgP);but Organization=OrgR U OrgS;thus Organization U OrgP=OrgR U (OrgS U OrgP)
				try:				OrgRUOrgP=TableUnion[(OrgR,OrgP)]									
				except KeyError:	continue
				try:				ResOrg=TableUnion[(OrgS,OrgRUOrgP)]
				except KeyError:	continue			
				TableUnion[(OrgP,Org)]=ResOrg
				#				TableUnion[(Org,OrgP)]=ResOrg				
				if ShowUnion:
				    print "%s U %s = (%s U %s) U %s = %s U (%s U %s) = %s U %s = %s"%(Org,OrgP,OrgS,OrgR,OrgP,OrgS,OrgR,OrgP,OrgS,OrgRUOrgP,ResOrg)
    			UnionByAUBUC2+=1    
    			return ResOrg
	# If the lattice was distributive here you could add that if A, B are such that A ^ B = Org then (A U OrgP) ^ (B U OrgP)=(A ^ B)U OrgP= Org U OrgP
	return None



	
	
	
	
	
	
	
	
def AddRelationAbove(Org,OrganizationsAboveNotCovering,OrganizationsMaybeCovering,lOrganization):
	Covering=set([])
	AddSetRelations(OrganizationsAboveNotCovering,set([Org]),(OBelowNotCovered,OAboveNotCovering))
	EliminateAnyRelations(OrganizationsAboveNotCovering,set([Org]))
	for OrgP in OrganizationsMaybeCovering:
		if NameSize[OrgP]==lOrganization+1:
			Covering|=set([OrgP])
			try:				OCovered[OrgP]|=set([Org])
			except KeyError:	OCovered[OrgP]=set([Org])
			try:				OCovering[Org]|=set([OrgP])
			except KeyError:	OCovering[Org]=set([OrgP])
		else:
			try:				OMaybeCovered[OrgP]|=set([Org])
			except KeyError:	OMaybeCovered[OrgP]=set([Org])
			try:				OMaybeCovering[Org]|=set([OrgP])
			except KeyError:	OMaybeCovering[Org]=set([OrgP])
	return Covering	
	
def EliminateTriangles(Org,LogOn=False):
	"sometimes as we add an organization in the middle we might shrtcircuit something. That is a relation that is above, was already connected with something that is below"
	if LogOn: LogThis("EliminateTriangles (%s)"%Org)
	OrgsNearBelow=OCovered[Org]|OMaybeCovered[Org]
	if LogOn: LogThis("OrgsNearBelow= %s"%OrgsNearBelow)
	OrgsNearAbove=OCovering[Org]|OMaybeCovering[Org]
	if LogOn: LogThis("OrgsNearAbove= %s"%OrgsNearAbove)
	AddSetRelations(OrgsNearAbove,OrgsNearBelow,(OBelowNotCovered,OAboveNotCovering))
	if LogOn: LogThis("AddRelations(%s,%s,3)"%(OrgsNearAbove,OrgsNearBelow))
	EliminateAnyRelations(OrgsNearAbove,OrgsNearBelow)
	if LogOn: LogThis("EliminateAnyRelations(%s,%s)"%(OrgsNearAbove,OrgsNearBelow))
	#	FindAllTriangles(OrgsNearAbove)	  ##########################################################################################
	return

def AddOrganizationUnionIntersection(Org):
	"""Here we are Trying to fill in the Two intersection and Union Table. We assume the Relation Table is complete. Still this is not enough to calculate all the union and intersection
	thus at times we find ourselves with new organizations coming up. It returns the set of organizations that could not solve"""
	OrgAb=OrganizationsAbove(Org)
	for OrgP in OrgAb:
		TableUnion[(Org, OrgP)]=OrgP
		#		TableUnion[(OrgP, Org)]=OrgP
		TableIntersection[(Org, OrgP)]=Org
		#		TableIntersection[(OrgP, Org)]=Org	
	OrgBel=OrganizationsUnder(Org)
	for OrgP in OrgBel:
		TableUnion[(Org, OrgP)]=Org
		#		TableUnion[(OrgP, Org)]=Org
		TableIntersection[(Org, OrgP)]=OrgP
		#		TableIntersection[(OrgP, Org)]=OrgP
	TableUnion[(Org,Org)]=Org
	TableIntersection[(Org,Org)]=Org
	ToBeSolvedByHandUnions=OSide[Org]|set([])
	ToBeSolvedByHandIntersections=OSide[Org]|set([])
	return ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections



def EliminateAnyRelations(OrgsA,OrgsB):
	#TODO what about the other relations? It assumes B>A
	for A in OrgsA:
		try:				OCovered[A]-=OrgsB
		except KeyError:	OCovered[A]=set([])
		try:				OMaybeCovered[A]-=OrgsB
		except KeyError:	OMaybeCovered[A]=set([])	
	for B in OrgsB:
		try:				OCovering[B]-=OrgsA
		except KeyError:	OCovering[B]=set([])
		try:				OMaybeCovering[B]-=OrgsA
		except KeyError:	OMaybeCovering[B]=set([])	


def OrganizationsUnder(Org):				return OCovered[Org]|OMaybeCovered[Org]|OBelowNotCovered[Org]
def OrganizationsAbove(Org):				return OCovering[Org]|OMaybeCovering[Org]|OAboveNotCovering[Org]
def OrganizationsBetween(TopOrg,BottomOrg):	return OrganizationsUnder(TopOrg)&OrganizationsAbove(BottomOrg)
def GetMaxOrganizationsFast(Os):
	"""Returns the Maximal subset of the set of Organizations. That is the set of organizations such that there is no organization lower than tham
	Cost=|Os|*|Os|*(|OAbove|+|OUnder|)
	"""
	Ol=list(Os)
	l=len(Ol)
	i=0
	j=1
	while i<l:
		while j<l:
			if Ol[j] in OrganizationsAbove(Ol[i]):
				l-=1
				Ol[i]=Ol[l]
				j=i+1
				break
			elif Ol[j] in OrganizationsUnder(Ol[i]):
				l-=1
				Ol[j]=Ol[l]
			else:
				j+=1
		else:
			i+=1
			j=i+1
	return set(Ol[:l])


def GetMinOrganizationsFast(Os):
	"""Returns the Minimal subset of the set of Organizations. That is the set of organizations such that there is no organization lower than tham
	Cost=|Os|*|Os|*(|OAbove|+|OUnder|)
	"""
	#thanks http://stackoverflow.com/questions/4117859/find-in-a-dynamic-pythonic-way-the-minimum-elements-in-a-partially-ordered-set/4117941#4117941
	Ol=list(Os)
	l=len(Ol)
	i=0
	j=1
	while i<l:
		while j<l:
			if Ol[j] in OrganizationsUnder(Ol[i]):
				l-=1
				Ol[i]=Ol[l]
				j=i+1
				break
			elif Ol[j] in OrganizationsAbove(Ol[i]):
				l-=1
				Ol[j]=Ol[l]
			else:
				j+=1
		else:
			i+=1
			j=i+1
	return set(Ol[:l])


def GetSetsRelationshipAny(Orgs,R,Context=OrganizationsStudied):
	Result=set([])
	for O in Orgs:
		Result|=R[O]
	return Result&Context


def AddRelationBelow(Org,OrganizationsBelowNotCovered,OrganizationsMaybeCovered,lOrganization):
	Covered=set([])
	AddSetRelations(set([Org]),OrganizationsBelowNotCovered,(OBelowNotCovered,OAboveNotCovering))
	EliminateAnyRelations(set([Org]),OrganizationsBelowNotCovered)
	for OrgP in OrganizationsMaybeCovered:
		if lOrganization==NameSize[OrgP]+1:
			Covered|=set([OrgP])
			try:				OCovered[Org]|=set([OrgP])
			except KeyError:	OCovered[Org]=set([OrgP])
			try:				OCovering[OrgP]|=set([Org])
			except KeyError:	OCovering[OrgP]=set([Org])
		else:
			try:				OMaybeCovered[Org]|=set([OrgP])
			except KeyError:	OMaybeCovered[Org]=set([OrgP])
			try:				OMaybeCovering[OrgP]|=set([Org])
			except KeyError:	OMaybeCovering[OrgP]=set([Org])
	return Covered	


def OrganizationsRelated2(TopOrg,BottomOrg):
	"""Cost:|1|"""
	OrganizationsOverTheTop=OrganizationsAbove(TopOrg)
	OrganizationsUnderTheBottom=OrganizationsUnder(BottomOrg)
	OrganizationsOnSideToBoth=OSide[TopOrg]&OSide[BottomOrg]
	organizationsbetween=OrganizationsBetween(TopOrg,BottomOrg)
	OrganizationsSideTop=OrganizationsAbove(BottomOrg)&OSide[TopOrg] 				#what is on the site the top, but above the bottom
	OrganizationsSideBottom=OSide[BottomOrg]&OrganizationsUnder(TopOrg)
	return (OrganizationsUnderTheBottom,OrganizationsSideBottom,OrganizationsOnSideToBoth,organizationsbetween,OrganizationsSideTop,OrganizationsOverTheTop)

def AddSetRelations(OrgsA,OrgsB,Relation):
	for A in OrgsA:
		try:				Relation[0][A]|=OrgsB
		except KeyError:	Relation[0][A]=OrgsB	
	for B in OrgsB:
		try:				Relation[1][B]|=OrgsA
		except KeyError:	Relation[1][B]=OrgsA

def AddSetRelation(A,B,Relation):
	try:				Relation[0][A]|=set([B])
	except KeyError:	Relation[0][A]=set([B])
	try:				Relation[1][B]|=set([A])
	except KeyError:	Relation[1][B]=set([A])


def AddRelationsByHand(Org,OrgsToCheck,lOrganization=None):
	if not OrgsToCheck: 								return
	if lOrganization==None:								lOrganization=NameSize[Org]
	OrganizationsBelowNotCovered=		set([])
	OrganizationsMaybeCovered=			set([])
	OrgsCovered=						set([])
	OrgsCovering=						set([])
	OrganizationsMaybeCovering=			set([])
	OrganizationsAboveNotCovering=		set([])
	OrgsAbove=							set([])
	OrgsBelow=							set([])
	OrgsUncomparable=					set([])
	for OrgP in OrgsToCheck:
		if 		IsOrg1inOrg2Simple( OrgP, Org):			OrgsBelow		|=set([OrgP])
		elif 	IsOrg1inOrg2Simple( Org, OrgP):			OrgsAbove		|=set([OrgP])
		else:											OrgsUncomparable|=set([OrgP])
	AddSetRelations(set([Org]),OrgsUncomparable,(OSide,OSide))
	if len(OrgsBelow):
		OrganizationsMaybeCovered=GetMaxOrganizationsFast(OrgsBelow)
		OrganizationsBelowNotCovered=OrgsBelow-OrganizationsMaybeCovered
		OrgsCovered=AddRelationBelow(Org,OrganizationsBelowNotCovered,OrganizationsMaybeCovered,lOrganization)
	if len(OrgsAbove):
		OrganizationsMaybeCovering=GetMinOrganizationsFast(OrgsAbove)
		OrganizationsAboveNotCovering=OrgsAbove-OrganizationsMaybeCovering
		OrgsCovering=AddRelationAbove(Org,OrganizationsAboveNotCovering,OrganizationsMaybeCovering,lOrganization)
	return (OrganizationsBelowNotCovered,OrganizationsMaybeCovered,OrgsCovered,OrgsUncomparable,OrgsCovering,OrganizationsMaybeCovering,OrganizationsAboveNotCovering)	




def AddOrganizationRelation(Org,TopOrg,BottomOrg): #This assumes that all the other relations before are there and also the two table union and intersection are complete
	lOrganization=len(NameOrganization.check(Org))	
	NameSize[Org]=lOrganization
	OCovered[Org]=set([])
	OMaybeCovered[Org]=set([])
	OCovering[Org]=set([])
	OMaybeCovering[Org]=set([])
	OAboveNotCovering[Org]=set([])
	OBelowNotCovered[Org]=set([])
	OSide[Org]=set([])
	LogThis("AddOrganizationRelation (%s, %s, %s)"%(Org,TopOrg,BottomOrg))
	(OrganizationsUnderTheBottom,OrganizationsSideBottom,OrganizationsOnSideToBoth,OrganizationsBetween,OrganizationsSideTop,OrganizationsOverTheTop)=OrganizationsRelated2(TopOrg,BottomOrg)
	AddSetRelations(set([Org]),OrganizationsUnderTheBottom,(OBelowNotCovered,OAboveNotCovering))
	AddSetRelations(set([Org]),OrganizationsOnSideToBoth,(OSide,OSide))
	AddSetRelations(OrganizationsOverTheTop,set([Org]),(OBelowNotCovered,OAboveNotCovering))
	LogThis("OrganizationsUnderTheBottom=%s"%OrganizationsUnderTheBottom)
	LogThis("OrganizationsSideBottom=%s"%OrganizationsSideBottom)
	LogThis("OrganizationsOnSideToBoth=%s"%OrganizationsOnSideToBoth)
	LogThis("OrganizationsBetween=%s"%OrganizationsBetween)
	LogThis("OrganizationsSideTop=%s"%OrganizationsSideTop)
	LogThis("OrganizationsOverTheTop=%s"%OrganizationsOverTheTop)
	OrgsAbove=set([TopOrg])
	OrgsBelow=set([BottomOrg])
	OrgsUncomparable=set([])
	for OrgP in OrganizationsBetween:
		if 		IsOrg1inOrg2Simple( OrgP, Org):			OrgsBelow		|=set([OrgP])
		elif 	IsOrg1inOrg2Simple( Org, OrgP):			OrgsAbove		|=set([OrgP])
		else:											OrgsUncomparable|=set([OrgP])
	AddSetRelations(set([Org]),OrgsUncomparable,(OSide,OSide))
	LogThis("OrgsAbove=%s"%OrgsAbove)
	LogThis("OrgsUncomparable=%s"%OrgsUncomparable)
	LogThis("BottomOrg=%s"%BottomOrg)
	OrganizationsMaybeCovering=GetMinOrganizationsFast(OrgsAbove)
	OrganizationsAboveNotCovering=OrgsAbove-OrganizationsMaybeCovering
	LogThis("OrganizationsMaybeCovering=%s"%OrganizationsMaybeCovering)
	LogThis("OrganizationsAboveNotCovering=%s"%OrganizationsAboveNotCovering)
	AddRelationAbove(Org,OrganizationsAboveNotCovering,OrganizationsMaybeCovering,lOrganization)
	OrganizationsMaybeCovered=GetMaxOrganizationsFast(OrgsBelow)
	OrganizationsBelowNotCovered=OrgsBelow-OrganizationsMaybeCovered
	AddRelationBelow(Org,OrganizationsBelowNotCovered,OrganizationsMaybeCovered,lOrganization)		
	OrgBelowComplicated=set([])
	for OrgP in OrganizationsSideBottom:
		#####################################TO CHECK IF IT IS USEFUL AT ALL###################################################################
		OrgUnionOrgPOrganization=TableUnion[(OrgP,BottomOrg)] #What if the table union is not complete? Shouldn't it just be a try?
		if OrgUnionOrgPOrganization in OrganizationsUnder(Org):
			AddSetRelation(Org,OrgP,(OBelowNotCovered,OAboveNotCovering))
		elif Org==OrgUnionOrgPOrganization:
			OrgBelowComplicated|=set([OrgP])
		else:											#if (OrgP U NameBottomOrg)  [not <] Org  =>  OrgP [not <] Org
			AddSetRelation(Org,OrgP,(OSide,OSide))
	if OrgBelowComplicated:
		SideBottomCovered=GetMaxOrganizationsFast(OrgBelowComplicated)
		SideBottomNotCovered=OrgBelowComplicated-SideBottomCovered
		for O in SideBottomCovered:
			if lOrganization==NameSize[O]+1:
				AddSetRelation(Org,O,(OCovering,OCovered))
			else:	
				AddSetRelation(Org,O,(OMaybeCovering,OMaybeCovered))
		for O in SideBottomNotCovered:
			AddSetRelation(Org,O,(OBelowNotCovered,OAboveNotCovering))
	OrgAboveComplicated=set([])
	for OrgP in OrganizationsSideTop:
		OrgIntersectionOrgPOrganization=TableIntersection[(OrgP,TopOrg)]
		if OrgIntersectionOrgPOrganization in OrganizationsAbove(Org):
			AddSetRelation(OrgP,Org,(OBelowNotCovered,OAboveNotCovering))
		elif OrgIntersectionOrgPOrganization==Org:
			OrgAboveComplicated|=set([OrgP])
		else:
			AddSetRelation(Org,OrgP,(OSide,OSide))
	if OrgAboveComplicated:
		SideTopCovering=GetMaxOrganizationsFast(OrgAboveComplicated)
		SideTopNotCovering=OrgAboveComplicated-SideTopCovering
		for O in SideTopCovering:
			if NameSize[O]==lOrganization+1:
				AddSetRelation(O,Org,(OCovering,OCovered))
			else:			
				AddSetRelation(O,Org,(OMaybeCovering,OMaybeCovered))
		for O in SideTopNotCovering:
			AddSetRelation(Org,O,(OBelowNotCovered,OAboveNotCovering))
	EliminateTriangles(Org)	
	global OrganizationsStudied
	OrganizationsStudied|=set([Org])


def AddOrganizationOfficially(Org,TopOrg,BottomOrg):
	global OrganizationsStudied
	AddOrganizationRelation(Org,TopOrg,BottomOrg)
	ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections=AddOrganizationUnionIntersection(Org)
	return ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections



















def AddOrganizationRelationFromUnion(Org,BottomOrg1,BottomOrg2): #This assumes that all the other relations before are there and also the two table union and intersection are complete
	try:
		lOrganization=len(NameOrganization.check(Org))	
	except IndexError:
		print "Problem with ",Org,": ",NameOrganization.check(Org)
		raise
	global OrganizationsStudied
	
	NameSize[Org]=lOrganization
	
	OCovered[Org]=set([])
	OMaybeCovered[Org]=set([])
	OCovering[Org]=set([])
	OMaybeCovering[Org]=set([])
	OAboveNotCovering[Org]=set([])
	OBelowNotCovered[Org]=set([])
	OSide[Org]=set([])	

	RelatedOrgs=set([Org])	

	OrgsAbove=(OMaybeCovering[BottomOrg1]|OAboveNotCovering[BottomOrg1])&(OMaybeCovering[BottomOrg2]|OAboveNotCovering[BottomOrg2])&(OrganizationsStudied-RelatedOrgs)
	OrganizationsMaybeCovering=GetMinOrganizationsFast(OrgsAbove)
	OrganizationsAboveNotCovering=OrgsAbove-OrganizationsMaybeCovering
	Covering=AddRelationAbove(Org,OrganizationsAboveNotCovering,OrganizationsMaybeCovering,lOrganization)
	RelatedOrgs|=OrgsAbove
	if len(Covering)>1:	
		TopOrg1,TopOrg2=list(Covering)[0:2]
	elif Covering==1 & len(OrganizationsMaybeCovering-Covering)>1:
		TopOrg1=list(Covering)[0]
		TopOrg2=list(OrganizationsMaybeCovering-Covering)[0]
	else:
		assert len(OrganizationsMaybeCovering)
		OrgsSideTop=GetSetsRelationshipAny(OrganizationsMaybeCovering,OSide,Context=OrganizationsStudied-RelatedOrgs)
		OrgsSide=GetSetsRelationshipAny((BottomOrg1,BottomOrg2),OSide,Context=OrgsSideTop)
		AddSetRelations(set([Org]),OrgsSide,(OSide,OSide))
		RelatedOrgs|=OrgsSide		
		OrgsBelowGen=(OrganizationsUnder(BottomOrg1)|OrganizationsUnder(BottomOrg2))&(OrganizationsStudied-RelatedOrgs)
		AddSetRelations(set([Org]),OrgsBelowGen,(OBelowNotCovered,OAboveNotCovering)) #I don't have to delete any previous relation, because Org is a new entry		
		RelatedOrgs|=OrgsBelowGen
		AddRelationsByHand(Org,OrganizationsStudied-RelatedOrgs,lOrganization)
		EliminateTriangles(Org,LogOn=True)	
		OrganizationsStudied|=set([Org])
		return
	OrgsSide=(OSide[TopOrg1]|OSide[TopOrg2])&(OSide[BottomOrg1]|OSide[BottomOrg2])&(OrganizationsStudied-RelatedOrgs)
	AddSetRelations(set([Org]),OrgsSide,(OSide,OSide))
	OrgsBelow=(OMaybeCovered[TopOrg1]|OBelowNotCovered[TopOrg1])&(OMaybeCovered[TopOrg2]|OBelowNotCovered[TopOrg2])&(OrganizationsStudied-RelatedOrgs)
	OrganizationsMaybeCovered=GetMaxOrganizationsFast(OrgsBelow)
	OrganizationsBelowNotCovered=OrgsBelow-OrganizationsMaybeCovered
	AddRelationBelow(Org,OrganizationsBelowNotCovered,OrganizationsMaybeCovered,lOrganization)
	EliminateTriangles(Org,LogOn=False)	
	OrganizationsStudied|=set([Org])	
	return

def AddOrganizationRelationFromIntersection(Org,TopOrg1,TopOrg2): #This assumes that all the other relations before are there and also the two table union and intersection are complete
	lOrganization=len(NameOrganization.check(Org))	
	global OrganizationsStudied
	NameSize[Org]=lOrganization
	OCovered[Org]=set([])
	OMaybeCovered[Org]=set([])
	OCovering[Org]=set([])
	OMaybeCovering[Org]=set([])
	OAboveNotCovering[Org]=set([])
	OBelowNotCovered[Org]=set([])
	OSide[Org]=set([])
	RelatedOrgs=set([Org]) #This will store all the organizations that has already been calculated
	OrgsBelow=(OMaybeCovered[TopOrg1]|OBelowNotCovered[TopOrg1])&(OMaybeCovered[TopOrg2]|OBelowNotCovered[TopOrg2])&(OrganizationsStudied-RelatedOrgs)
	OrganizationsMaybeCovered=GetMaxOrganizationsFast(OrgsBelow)
	OrganizationsBelowNotCovered=OrgsBelow-OrganizationsMaybeCovered
	Covered=AddRelationBelow(Org,OrganizationsBelowNotCovered,OrganizationsMaybeCovered,lOrganization)
	RelatedOrgs|=OrgsBelow
	if len(Covered)>1:	
		BottomOrg1,BottomOrg2=list(Covered)[0:2]
	elif Covered==1 & len(OrganizationsMaybeCovered-Covered)>1:
		BottomOrg1=list(Covered)[0]
		BottomOrg2=list(OrganizationsMaybeCovered-Covered)[0]
	else:		
		assert len(OrganizationsMaybeCovered)
		OrgsSideBottom=GetSetsRelationshipAny(OrganizationsMaybeCovered,OSide,Context=OrganizationsStudied-RelatedOrgs)
		OrgsSide=GetSetsRelationshipAny((TopOrg1,TopOrg2),OSide,Context=OrgsSideBottom)
		AddSetRelations(set([Org]),OrgsSide,(OSide,OSide)) #I don't have to delete any previous relation, because Org is a new entry		
		RelatedOrgs|=OrgsSide
		OrgsAboveGen=(OrganizationsAbove(TopOrg1)|OrganizationsAbove(TopOrg2))&(OrganizationsStudied-RelatedOrgs)
		AddSetRelations(OrgsAboveGen,set([Org]),(OBelowNotCovered,OAboveNotCovering)) #I don't have to delete any previous relation, because Org is a new entry		
		RelatedOrgs|=OrgsAboveGen
		AddRelationsByHand(Org,OrganizationsStudied-RelatedOrgs,lOrganization)
		EliminateTriangles(Org)	
		OrganizationsStudied|=set([Org])
		return
	OrgsSide=(OSide[TopOrg1]|OSide[TopOrg2])&(OSide[BottomOrg1]|OSide[BottomOrg2])&(OrganizationsStudied-RelatedOrgs)
	AddSetRelations(set([Org]),OrgsSide,(OSide,OSide))
	OrgsAbove=(OMaybeCovering[BottomOrg1]|OAboveNotCovering[BottomOrg1])&(OMaybeCovering[BottomOrg2]|OAboveNotCovering[BottomOrg2])&(OrganizationsStudied-RelatedOrgs)	
	OrganizationsMaybeCovering=GetMinOrganizationsFast(OrgsAbove)
	OrganizationsAboveNotCovering=OrgsAbove-OrganizationsMaybeCovering
	AddRelationAbove(Org,OrganizationsAboveNotCovering,OrganizationsMaybeCovering,lOrganization)
	EliminateTriangles(Org)
	OrganizationsStudied|=set([Org])
	return




def AddUnionToTable(OrgL, OrgR, OrgUnion):
	global UnionByAUBBuilding
	if OrgL==OrgUnion: return
	if OrgR==OrgUnion: return
	OrgsBelowClosure=OrganizationsUnder(OrgUnion)
	OrgsAboveL=OrganizationsAbove(OrgL)
	OrgsAboveR=OrganizationsAbove(OrgR)
	OrgsInterestingAboveL=(OrgsAboveL&OrgsBelowClosure)
	OrgsInterestingAboveR=(OrgsAboveR&OrgsBelowClosure)
	OrgsInterestingAboveL|=set([OrgL])
	OrgsInterestingAboveR|=set([OrgR])
	for OL in OrgsInterestingAboveL:
		for OR in OrgsInterestingAboveR:
			try:
				OrgLUOrgR=TableUnion[(OL,OR)]
				if (OrgUnion!=OrgLUOrgR):
					print "ERROR 4327842378924234893247983249843289"
					print "OrgL=",OrgL, NameSize[OrgL]
					print "OrgR=",OrgR, NameSize[OrgR]
					print "OrgUnion=",OrgUnion, NameSize[OrgUnion]
					print "OL=",OL, NameSize[OL]
					print "OR=",OR, NameSize[OR]
					print "OrgLUOrgR=",OrgLUOrgR, NameSize[OrgLUOrgR]
					closedatabases()
					sys.exit()									
			except KeyError:
				TableUnion[(OR,OL)]=OrgUnion   
				#				TableUnion[(OL,OR)]=OrgUnion
				UnionByAUBBuilding+=1

def AddIntersectionToTable(OrgL, OrgR, OrgIntersection):
	global IntersectionByAIBBuilding
	OrgsAboveIntersection=OrganizationsAbove(OrgIntersection)
	OrgsBelowL=OrganizationsUnder(OrgL)
	OrgsBelowR=OrganizationsUnder(OrgR)
	OrgsInterestingBelowL=(OrgsBelowL&OrgsAboveIntersection)
	OrgsInterestingBelowR=(OrgsBelowR&OrgsAboveIntersection)
	OrgsInterestingBelowL|=set([OrgL])
	OrgsInterestingBelowR|=set([OrgR])
	for OL in OrgsInterestingBelowL:
		for OR in OrgsInterestingBelowR:
			try:
				OrgLIOrgR=TableIntersection[(OL,OR)]
				if (OrgIntersection!=OrgLIOrgR):
					print "ERROR 432784sanjd47983249843289"
					sys.exit()									
			except KeyError:
				TableIntersection[(OR,OL)]=OrgIntersection   
				#				TableIntersection[(OL,OR)]=OrgIntersection
				IntersectionByAIBBuilding+=1
	







def AddOrganizationOfficiallyFromUnion(Org,BottomOrg1,BottomOrg2):
	global OrganizationsStudied
	AddOrganizationRelationFromUnion(Org,BottomOrg1,BottomOrg2)
	ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections=AddOrganizationUnionIntersection(Org)
	return ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections

def AddOrganizationOfficiallyFromIntersection(Org,TopOrg1,TopOrg2):
	global OrganizationsStudied
	AddOrganizationRelationFromIntersection(Org,TopOrg1,TopOrg2)
	ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections=AddOrganizationUnionIntersection(Org)
	return ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections



def BuildSetsBigger(SetsSmaller, SizeSets):
	"""This function is a little jewel. The idea is that we need to build sets of size n starting from sets of size n-1. 
	But we are allowed to do so only if there are all the subsets of size n-1 contained in the set of size n that we are about to build"""
	NumberSets	=	len(SetsSmaller)
	if 		NumberSets	==	0:				return set([])
	n =	SizeSets+1 													#equivalent to say: SizeSets = n-1
	if 		NumberSets	<	n:				return set([]) 			#to build a set of size n we need at least n sets of size n-1
	PossibleUnions=findsubsets(SetsSmaller,2)
	TimesFound=dict() 												#we count how many times each set of size n has been found
	for a in PossibleUnions:
		l=list(a)
		if SizeSets==1:
			NewSet=frozenset([a[0],a[1]])
			try:							TimesFound[NewSet]+=1
			except KeyError:				TimesFound[NewSet] =1
		elif len(l[0] & l[1])==SizeSets-1:
			NewSet = a[0] | a[1]
			try: 							TimesFound[NewSet]+=1
			except KeyError:				TimesFound[NewSet] =1
	Unions=set([])
	NAim=(n*(n-1))/2
	for t in TimesFound.keys():
		if TimesFound[t]==NAim:
			Unions	|=		set([t])
	return Unions


def FindNewOrganizationsAboveByY(BottomOrg,molecules,TopOrg,OrganizationsKnown,SizeSets=2,MaxNumberOrg=None):
	"""Given an organizations BottomOrg, 
	finds all the organizations generated by bottomorg adding SizeSets molecules at a time.
	and each time it finds an organisation it expands the sub-lattice
	"""
	
	StartTimeX=time.clock()

	BottomOrgDetails=NameOrganization.check(BottomOrg)
	
	print "I am testing groups of size ",SizeSets,"above the org ",BottomOrg, " (",BottomOrgDetails,")"
	
	PreviousSize=SizeSets-1
	try:											
		#ToBeCheckedMolecules=DownwardMolecules[BottomOrg][PreviousSize]
		ToBeCheckedSets=DownwardMolecules[BottomOrg][PreviousSize]
	except KeyError:								return set([])
	#MoleculesTodo=len(ToBeCheckedSets)
	SetsTodo=len(ToBeCheckedSets)
	NewOrganizations=set([])
	#ThisLayerOrganizations=set([])
	print "Looking for the organisations above ",BottomOrg,"; Sets to be checked: ",SetsTodo,": ",ToBeCheckedSets
	print "  UpwardSets:",len(  UpwardMolecules[BottomOrg][PreviousSize]),":",  UpwardMolecules[BottomOrg][PreviousSize]
	print "SidewardSets:",len(SidewardMolecules[BottomOrg][PreviousSize]),":",SidewardMolecules[BottomOrg][PreviousSize]
	print "DownwardSets:",len(DownwardMolecules[BottomOrg][PreviousSize]),":",DownwardMolecules[BottomOrg][PreviousSize]
	#missing the DiagonalSets (the ones where only part of the molecules are inside)

	if SizeSets not in DownwardMolecules[BottomOrg].keys():         DownwardMolecules[BottomOrg][SizeSets]=set([])
	if SizeSets not in   UpwardMolecules[BottomOrg].keys():           UpwardMolecules[BottomOrg][SizeSets]=set([])
	if SizeSets not in SidewardMolecules[BottomOrg].keys():         SidewardMolecules[BottomOrg][SizeSets]=set([])

	print "  UpwardMolecules:[",SizeSets,"]",len(  UpwardMolecules[BottomOrg][SizeSets]),":",  UpwardMolecules[BottomOrg][SizeSets]
	print "SidewardMolecules:[",SizeSets,"]",len(SidewardMolecules[BottomOrg][SizeSets]),":",SidewardMolecules[BottomOrg][SizeSets]
	print "DownwardMolecules:[",SizeSets,"]",len(DownwardMolecules[BottomOrg][SizeSets]),":",DownwardMolecules[BottomOrg][SizeSets]


	#	DeadEnd=set()
	#	DeadEnd=SidewardMolecules[BottomOrg][SizeSets-1]|UpwardMolecules[BottomOrg][SizeSets-1]

	OriginalOrganizationsFound=0

	print "ToBeCheckedSets=",ToBeCheckedSets
	
	SetsToTest=BuildSetsBigger(ToBeCheckedSets,SizeSets-1)
	for a in SetsToTest:
		SetsTodo-=1
		ThisLayerOrganizations=set([])	
		SetA=set(a)	
		print "Looking for the organisations above ",BottomOrg,"(",BottomOrgDetails,")","; Set to be checked: ",SetA,
		S=SetA|BottomOrgDetails
											#TODO here I should check if the set is directly an organisation known, before calculating it
		OrgS=Organization(S)
		if OrgS==BottomOrgDetails:
			DownwardMolecules[BottomOrg][SizeSets].add(a)
			print
			continue
		else:
			NameS=OrganizationName.check(OrgS)						#			ExplosiveSets|=set([frozenset([a])])
			print " Org: %s := G(%sU{%s})=%s"	%(NameS,BottomOrg,a,OrgS)
			
			if NameS not in NewOrganizations|OrganizationsKnown:
				#				print "New Org: ", NameS, ":= G(", BottomOrg, "U{",a,"}) = "	
				
				print "New Org: %s := G(%sU{%s})=%s"	%(NameS,BottomOrg,a,OrgS)
				OriginalOrganizationsFound+=1
				NewOrganizations|=set([NameS])
				#ThisLayerOrganizations|=set([NameS])
				print "adding", NameS
				sys.stdout.flush()
				
				ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections=AddOrganizationOfficially(NameS,TopOrg,BottomOrg)
				print NameS, "added"
				sys.stdout.flush()
				
				if ToBeSolvedByHandUnions:					ToBeSolvedByHandUnionsDict[NameS]=ToBeSolvedByHandUnions
				if ToBeSolvedByHandIntersections:			ToBeSolvedByHandIntersectionsDict[NameS]=ToBeSolvedByHandIntersections
				print DescribeSituation(),
				sys.stdout.flush()
				
				#we found an organisation, so now we complete the lattice
				while ToBeSolvedByHandUnionsDict or ToBeSolvedByHandIntersectionsDict:
					
					keyslength=[(k,NameSize[k]) for k in ToBeSolvedByHandUnionsDict.keys()]
					keyslength.sort(sortpairs, reverse=False)
					for k,s in keyslength:
						#print "k=",k,
						#print "s=",s
						NewOrganizations|=SolveByHandUnion(k)
						
					keyslength=[(k,NameSize[k]) for k in ToBeSolvedByHandIntersectionsDict.keys()]
					keyslength.sort(sortpairs, reverse=True)
					for k,s in keyslength:
						NewOrganizations|=SolveByHandIntersection(k)
										
				print
				print DescribeSituation(),
				print " sublattice",
				#				buf=TestAllRelations();
				#				if buf:			print buf
				#				else: print "all relations are fine."
				if MaxNumberOrg:
					if len(NewOrganizations)>MaxNumberOrg:
						return NewOrganizations
						
						#closedatabases()
						#sys.exit(0)

			if a in OrgS:
 				UpwardMolecules[BottomOrg][SizeSets].add(a)
				for o in OrganizationsBetween(NameS,BottomOrg):
					try:                UpwardMolecules[o][SizeSets].add(a)
					except KeyError:    UpwardMolecules[o][SizeSets]=set([a])
			else:
				SidewardMolecules[BottomOrg][SizeSets].add(a)
				
				if NameS not in DownwardMolecules.keys():           
					DownwardMolecules[NameS]=dict()
					DownwardMolecules[NameS][SizeSets]=set([a])
				else:
					try:				DownwardMolecules[NameS][SizeSets].add(a)
					except KeyError:	DownwardMolecules[NameS][SizeSets]=set([a])
				if NameS not in   UpwardMolecules.keys():           UpwardMolecules[NameS]  =dict()
				if NameS not in SidewardMolecules.keys():           SidewardMolecules[NameS]=dict()
				
				for o in OrganizationsBetween(NameS,BottomOrg):
					if o not in   UpwardMolecules.keys():           
						UpwardMolecules[o]  		=dict()
						UpwardMolecules[o][SizeSets]=set([a])
					else:
						try:				UpwardMolecules[o][SizeSets].add(a)
						except KeyError:	UpwardMolecules[o][SizeSets]=set([a])
					if o not in SidewardMolecules.keys():           SidewardMolecules[o]=dict()
					if o not in DownwardMolecules.keys():           DownwardMolecules[o]=dict()
					
	print DescribeSituation()
	sys.stdout.flush()
	return NewOrganizations




def FindNewOrganizationsAbove(BottomOrg,molecules,TopOrg,OrganizationsKnown,MaxNumberOrg=None):
	"""Given an organizations BottomOrg, 
	finds all the organizations generated by bottomorg adding one molecule at a time.
	and each time it finds an organisation it expands the sub-lattice
	"""
	BottomOrgDetails=NameOrganization.check(BottomOrg)
	ToBeCheckedMolecules=molecules-BottomOrgDetails
	MoleculesTodo=len(ToBeCheckedMolecules)
	NewOrganizations=set([])
	#ThisLayerOrganizations=set([])
	print "Looking for the organisations above ",BottomOrg,"; Molecules to be checked: ",MoleculesTodo,": ",ToBeCheckedMolecules
	#if BottomOrg not in DownwardMolecules.keys():           DownwardMolecules[BottomOrg]=set([])
	if BottomOrg not in DownwardMolecules.keys():           DownwardMolecules[BottomOrg]=dict()
	if BottomOrg not in   UpwardMolecules.keys():           UpwardMolecules[BottomOrg]  =dict()
	if BottomOrg not in SidewardMolecules.keys():           SidewardMolecules[BottomOrg]=dict()

	if 1 not in DownwardMolecules[BottomOrg].keys():         DownwardMolecules[BottomOrg][1]=set([])
	if 1 not in UpwardMolecules[BottomOrg].keys():           UpwardMolecules[BottomOrg][1]  =set([])
	if 1 not in SidewardMolecules[BottomOrg].keys():         SidewardMolecules[BottomOrg][1]=set([])
	
	
	print "  UpwardMolecules [1]:",len(  UpwardMolecules[BottomOrg][1]),":",  UpwardMolecules[BottomOrg][1]
	print "SidewardMolecules [1]:",len(SidewardMolecules[BottomOrg][1]),":",SidewardMolecules[BottomOrg][1]
	print "DownwardMolecules [1]:",len(DownwardMolecules[BottomOrg][1]),":",DownwardMolecules[BottomOrg][1]

	OriginalOrganizationsFound=0

	MolsCovered=DownwardMolecules[BottomOrg][1]|UpwardMolecules[BottomOrg][1]|SidewardMolecules[BottomOrg][1]
	ToBeCheckedMolecules=ToBeCheckedMolecules-MolsCovered
	for a in ToBeCheckedMolecules:
		MoleculesTodo-=1
		ThisLayerOrganizations=set([])		
		print "Looking for the organisations above ",BottomOrg,"; Molecule to be checked: ",a,
		S=set([a])|BottomOrgDetails
											#TODO here I should check if the set is directly an organisation known, before calculating it
		OrgS=Organization(S)
		if OrgS==BottomOrgDetails:
			DownwardMolecules[BottomOrg][1].add(a)
			print
			continue
		else:
			NameS=OrganizationName.check(OrgS)						#			ExplosiveSets|=set([frozenset([a])])
			print " Org: %s := G(%sU{%s})=%s"	%(NameS,BottomOrg,a,OrgS)
			if NameS not in NewOrganizations|OrganizationsKnown:
				#				print "New Org: ", NameS, ":= G(", BottomOrg, "U{",a,"}) = "	
				
				print "New Org: %s := G(%sU{%s})=%s"	%(NameS,BottomOrg,a,OrgS)
				OriginalOrganizationsFound+=1
				NewOrganizations|=set([NameS])
				#ThisLayerOrganizations|=set([NameS])
				print "adding", NameS
				sys.stdout.flush()
				
				ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections=AddOrganizationOfficially(NameS,TopOrg,BottomOrg)
				print NameS, "added"
				sys.stdout.flush()
				
				if ToBeSolvedByHandUnions:					ToBeSolvedByHandUnionsDict[NameS]=ToBeSolvedByHandUnions
				if ToBeSolvedByHandIntersections:			ToBeSolvedByHandIntersectionsDict[NameS]=ToBeSolvedByHandIntersections
				print DescribeSituation(),
				sys.stdout.flush()
				
				#we found an organisation, so now we complete the lattice
				while ToBeSolvedByHandUnionsDict or ToBeSolvedByHandIntersectionsDict:
					
					keyslength=[(k,NameSize[k]) for k in ToBeSolvedByHandUnionsDict.keys()]
					keyslength.sort(sortpairs, reverse=False)
					for k,s in keyslength:
						#print "k=",k,
						#print "s=",s
						NewOrganizations|=SolveByHandUnion(k)
						
					keyslength=[(k,NameSize[k]) for k in ToBeSolvedByHandIntersectionsDict.keys()]
					keyslength.sort(sortpairs, reverse=True)
					for k,s in keyslength:
						NewOrganizations|=SolveByHandIntersection(k)
										
				print
				print DescribeSituation(),
				print " sublattice",
				#				buf=TestAllRelations();
				#				if buf:			print buf
				#				else: print "all relations are fine."
				if MaxNumberOrg:
					if len(NewOrganizations)>MaxNumberOrg:
						return NewOrganizations
						
						#closedatabases()
						#sys.exit(0)
			if a in OrgS:
 				UpwardMolecules[BottomOrg][1].add(a)
				for o in OrganizationsBetween(NameS,BottomOrg):
					if o not in   UpwardMolecules.keys():           
						UpwardMolecules[o]   =dict()
						UpwardMolecules[o][1]=set([a])
					else:
						try:                UpwardMolecules[o][1].add(a)
						except KeyError:    UpwardMolecules[o][1]=set([a])
			else:
				SidewardMolecules[BottomOrg][1].add(a)
				if NameS not in   DownwardMolecules.keys():
					DownwardMolecules[NameS]   =dict()
					DownwardMolecules[NameS][1]=set([a])				
				else:
					try:				DownwardMolecules[NameS][1].add(a)
					except KeyError:	DownwardMolecules[NameS][1]=set([a])
				
				if NameS not in   UpwardMolecules.keys():           UpwardMolecules[NameS]  =dict()
				if NameS not in SidewardMolecules.keys():           SidewardMolecules[NameS]=dict()
				
				for o in OrganizationsBetween(NameS,BottomOrg):
					if o not in   UpwardMolecules.keys():
						UpwardMolecules[o]   =dict()
						UpwardMolecules[o][1]=set([a])
					else:
						try:				UpwardMolecules[o][1].add(a)
						except KeyError:	UpwardMolecules[o][1]=set([a])
			#for k in   UpwardMolecules.keys():
			#	if k != BottomOrg:					print "UpwardMolecules[%s]=%s"%(k,UpwardMolecules[k])
			#for k in   SidewardMolecules.keys():
			#	if k != BottomOrg:					print "SidewardMolecules[%s]=%s"%(k,SidewardMolecules[k])
			#for k in   DownwardMolecules.keys():
			#	if k != BottomOrg:					print "DownwardMolecules[%s]=%s"%(k,DownwardMolecules[k])	
	print DescribeSituation()
	sys.stdout.flush()
	return NewOrganizations


def MapLatticeFromBelow(NameTopOrg, NameBottomOrg,MaxNumberOrg=None):
	molecules=NameOrganization.check(NameTopOrg)
	assert IsOrg1inOrg2Simple( NameBottomOrg,NameTopOrg)
	if NameTopOrg not in NameSize:		
		NameSize[NameTopOrg]=		len(NameOrganization.check(NameTopOrg))
		OCovered[NameTopOrg]=set([])
		OMaybeCovered[NameTopOrg]=set([])
		OCovering[NameTopOrg]=set([])
		OMaybeCovering[NameTopOrg]=set([])
		OAboveNotCovering[NameTopOrg]=set([])
		OBelowNotCovered[NameTopOrg]=set([])
		OSide[NameTopOrg]=set([])
		
	if NameBottomOrg not in NameSize:	
		NameSize[NameBottomOrg]=	len(NameOrganization.check(NameBottomOrg))
		OCovered[NameBottomOrg]=set([])
		OMaybeCovered[NameBottomOrg]=set([])
		OCovering[NameBottomOrg]=set([])
		OMaybeCovering[NameBottomOrg]=set([])
		OAboveNotCovering[NameBottomOrg]=set([])
		OBelowNotCovered[NameBottomOrg]=set([])
		OSide[NameBottomOrg]=set([])
		
	global OrganizationsStudied
	OrganizationsStudied|=set([NameTopOrg])
	OrganizationsStudied|=set([NameBottomOrg])
		
	AddSetRelation(NameTopOrg,NameBottomOrg,(OMaybeCovered,OMaybeCovering))
	TableUnion[(NameTopOrg,NameBottomOrg)]=NameTopOrg
	TableUnion[(NameTopOrg,NameTopOrg)]=NameTopOrg
	TableUnion[(NameBottomOrg,NameBottomOrg)]=NameBottomOrg
	TableIntersection[(NameTopOrg,NameBottomOrg)]=NameBottomOrg
	TableIntersection[(NameBottomOrg,NameBottomOrg)]=NameBottomOrg
	TableIntersection[(NameTopOrg,NameTopOrg)]=NameTopOrg

	OrganizationsKnown=set([]);
	OrganizationsKnown|=set([NameTopOrg])
	OrganizationsKnown|=set([NameBottomOrg])
	
	OrganisationToExplore=2
	global starttime
	starttime=time.clock()
	

	OrganizationsUsedAsBase=set([]);
	DepthOfExplorationByOrganisation=dict() #we assign to each organisation a value indicating what should be done next. 0= nothing, 1 explore 1 molecule at a time, 2, explore 2 molecules at a time, ...
	DepthOfExplorationByOrganisation[1]=0 #Organisation 1 is the top organisation, and makes no sense to go up.
	DepthOfExplorationByOrganisation[2]=1
	while sum(DepthOfExplorationByOrganisation.values()):
		for k in DepthOfExplorationByOrganisation.keys():
			NextAction=DepthOfExplorationByOrganisation[k]
			if NextAction==1:
				LastOrgsFound=FindNewOrganizationsAbove(k,molecules,NameTopOrg,OrganizationsKnown,MaxNumberOrg)
				for NewOrg in LastOrgsFound:
					try:
						DepthOfExplorationByOrganisation[NewOrg]
					except KeyError:
						DepthOfExplorationByOrganisation[NewOrg]=1
						OrganizationsKnown.add(NewOrg)
						if MaxNumberOrg:
							if len(OrganizationsKnown)>MaxNumberOrg:
								closedatabases()
								sys.exit(0)
				if len(DownwardMolecules[k][NextAction])>NextAction: 	#you need at least n sets of size n-1 to make a set of size n
																		DepthOfExplorationByOrganisation[k]=2
				else:													DepthOfExplorationByOrganisation[k]=0
			elif NextAction>1:
				LastOrgsFound=FindNewOrganizationsAboveByY(k,molecules,NameTopOrg,OrganizationsKnown,NextAction,MaxNumberOrg)
				for NewOrg in LastOrgsFound:
					try:	DepthOfExplorationByOrganisation[NewOrg]
					except KeyError:
							DepthOfExplorationByOrganisation[NewOrg]=1
							OrganizationsKnown.add(NewOrg)
							if MaxNumberOrg:
								if len(OrganizationsKnown)>MaxNumberOrg:
									closedatabases()
									sys.exit(0)
				if len(DownwardMolecules[k][NextAction])>NextAction:	#you need at least n sets of size n-1 to make a set of size n
																		DepthOfExplorationByOrganisation[k]=NextAction+1
				else:													DepthOfExplorationByOrganisation[k]=0
		
	print "Organizations Known=",
	for OrgThis in DepthOfExplorationByOrganisation.keys():
		print "Organisation ",OrgThis,"=",NameOrganization.check(OrgThis)
	
	

