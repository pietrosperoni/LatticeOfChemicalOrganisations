# -*- coding: utf-8 -*-
import os
#import sys
import gc
import random
#import math
#import copy


import time

import cPickle
import pickle

#import shelve

#import itertools


### code is poetry ###

###########################################
######## Constants ########################
###########################################

os.environ['PATH'] += ":"+"/usr/local/bin"
from subprocess import check_call

Org2FillColor={}

################################################################################

def linearapproximation(minOvalue,maxOvalue,value,minIvalue=0,maxIvalue=1):
	value=(float)(value-minIvalue)/maxIvalue
	delta=maxOvalue-minOvalue
	return (value*delta)+minOvalue

	
def next_string(s):
	"returns the next alphabetic string, only works with lowercase"
	strip_zs = s.rstrip('z')
	if strip_zs:
		return strip_zs[:-1] + chr(ord(strip_zs[-1]) + 1) + 'a' * (len(s) - len(strip_zs))
	else:
		return 'a' * (len(s) + 1)



def sumsequence(sequence,offset=0):
	result=list(sequence[:])
	result[0]+=offset
	for t in xrange(1,len(result)):
		result[t]+=result[t-1]
	return result


def find_keys(dic, val):
	"""return the key of dictionary dic given the value"""
	return [k for k, v in symbol_dic.iteritems() if v == val]


def ValMolecules1(sequence): return sequence
def ValMolecules2(sequence): return sumsequence(sequence)
def LayoutDrawType():
	while 1:
		yield ValMolecules1
		yield ValMolecules2
		
def findconvered(s,S):
	coveringrelations=[]
	for t in S:
		if s == t:						continue
		if t < s:
			for r in S:
				if t < r < s:
					break
			else:
				coveringrelations.append((s,t))
	return coveringrelations



def FindCovers(S):
	coverrelations=[]
	for s in S:
		a=findconvered(s,S)
		for rel in a:
			coverrelations.append(rel)
	return coverrelations
		


def WriteLattice(S,directory,name,sublattice=set([]), Org2Color={}):
	relations=FindCovers(S)
	buf="digraph G {%s%s"%(os.linesep,os.linesep)
	OrgBySize={}
	for s in S:
		try:
			OrgBySize[len(s)]=OrgBySize[len(s)]+[s]
		except KeyError:
			OrgBySize[len(s)]=[s]
	
	layers=[]
	r0=""
	
	for s in OrgBySize.values():
		if len(s)<1:
			continue
		r0='{rank=same; '
		#		 r0+="%s"%OrganizationName.check(s)
		#		 r0+="\\n"
		#		 print s
		
		for o2 in s:
			o=list(o2)
			o.sort()
			
			r0+=' "'
			t=0
			l=len(o)
			width=math.ceil(l**.5)		   
			for molecule in o:
				r0+="%s"%molecule
				t+=1
				if t<l:
					if (t%width):
						r0+=" "
					else:
						r0+="\\n"
			if not l:
				r0+='∅'
				
			r0+='" '
		r0+=' "'
		r0+='%s'%len(s[0])
		layers.append(len(s[0]))
		r0+='" '
		
		r0+='};%s'%os.linesep
		buf+=r0
		
	layers.sort()
	layers.reverse()
	
	for lay in range(0,len(layers)-1):
		r0+='"%s" -> "%s" [color=white];%s'%(layers[lay],layers[lay+1],os.linesep)
	buf+=r0
	
	for lay in range(0,len(layers)):
		r0+='"%s" [color=white];%s'%(layers[lay],os.linesep)
	buf+=r0
	
	for r in relations:
		r0=r1='"'
		if r[0]:
			rlist=list(r[0])
			rlist.sort()
			
			t=0
			l=len(r[0])
			width=math.ceil(l**.5)		   
			#			 for molecule in r[0]:
			for molecule in rlist:
				
				r0+="%s"%molecule
				t+=1
				if t<l:
					if (t%width):
						r0+=" "
					else:
						r0+="\\n"
					
		else:
			r0+="∅"
		if r[1]:
			rlist=list(r[1])
			rlist.sort()
			
			t=0
			l=len(r[1])
			width=math.ceil(l**.5) 
			for molecule in rlist:
			#			 for molecule in r[1]:
				
				r1+="%s"%molecule
				t+=1
				if t<l:
					if (t%width):
						r1+=" "
					else:
						r1+="\\n"
		else:
			r1+="∅"
		r0+='"'
		r1+='"'
		
		##		  if set(r)<set(sublattice):
		##			  r2=" [color=red]"
		##		  else:
		##			  r2=""
		
		print "r[0]=",r[0]
		if r[0] in Org2Color:
			r2=' [color="%s"]'%Org2Color[r[0]]
		elif r[0] in Org2FillColor:
			r2=' [color="%s"]'%Org2FillColor[r[0]]
		else:
			r2=''
			
			
		buf+=r0+" -> "+r1+r2+"; %s"%os.linesep
		
		
	for o2 in sublattice:#colors the sublattice red
		o=list(o2)
		o.sort()
		r0='"'
		t=0
		l=len(o)
		width=math.ceil(l**.5)		   
		for molecule in o:
			r0+="%s"%molecule
			t+=1
			if t<l:
				if (t%width):
					r0+=" "
				else:
					r0+="\\n"
		if not l:
			r0+='∅'			   
		r0+='" [fillcolor=red, style=filled];%s'%os.linesep
		buf+=r0
		
		
	for o2 in Org2FillColor.keys():
		o=list(o2)
		o.sort()
		r0='"'
		t=0
		l=len(o)
		width=math.ceil(l**.5)		   
		for molecule in o:
			r0+="%s"%molecule
			t+=1
			if t<l:
				if (t%width):
					r0+=" "
				else:
					r0+="\\n"
		if not l:
			r0+='∅'
		if o2 in Org2Color:
			r0+='" [fillcolor="%s", peripheries=6, color="%s",style=filled];%s'%(Org2FillColor[o2],Org2Color[o2],os.linesep)
		else:
			r0+='" [fillcolor="%s", style=filled];%s'%(Org2FillColor[o2],os.linesep)			
		buf+=r0
	
	buf+="}%s%s"%(os.linesep,os.linesep)
	f=open("%s%s%s.dot"%(directory,os.sep,name),'w')
	f.write(buf)
	f.close()

def WriteLatticesKnownInGeneral():
		AllOrg=GetAllKnownOrg()
		WriteLattice(AllOrg,'.',"L_AllOrganizationsGeneral")
		print "writing the graph"
		print check_call(["dot", "-ogeneral.png", "-Tpng", "./L_AllOrganizationsGeneral.dot"])
		print "graph drawn"

	
################################################################################

def LevDistance(a,b):
	"""Calculates the Levenshtein distance between a and b.
	From http://hetland.org/coding/python/levenshtein.py"""
	n, m = len(a), len(b)
	if n > m:
		# Make sure n <= m, to use O(min(n,m)) space
		a,b = b,a
		n,m = m,n
		
	current = range(n+1)
	for i in range(1,m+1):
		previous, current = current, [i]+[0]*n
		for j in range(1,n+1):
			add, delete = previous[j]+1, current[j-1]+1
			change = previous[j-1]
			if a[j-1] != b[i-1]:
				change = change + 1
			current[j] = min(add, delete, change)
			
	return current[n]


def updatearchive(step,b):
	print "I am updating the archive"
	AllStates["%10s"%step]=list(b.molecules)
	AllRndSeeds["%10s"%step]=random.getstate()


profiling=0
timelimit=0

def testtime(where,exceptionname):
	if profiling:
		if time.time()>timelimit:
			print where
			raise exceptionname



###########################################
######## Molecule Library #################
###########################################

def GenerateBasicOrganizations():
	#	 OrgLibrary.check(ReturnAllMolecules())
	OrgLibrary.check(set([])) 





def syncdatabases():
	global moleculelibrary, OrgLibrary, reactionDB, mutationsDB, LatticeLibrary,TimeBiggerOrganizationName,IntersectionRelOrg
	moleculelibrary.DB.sync()
	OrgLibrary.DB.sync()
	reactionDB.DB.sync()
	mutationsDB.DB.sync()
	LatticeLibrary.DB.sync()
	GeneratorsLibrary.DB.sync()
	OrganizationName.DB.sync()
	NameOrganization.DB.sync()
	TimeBiggerOrganizationName.DB.sync()
	IntersectionRelOrg.DB.sync()		
	global AllStates, AllRndSeeds
	print "I am sync the archives"
	AllStates.sync()
	AllRndSeeds.sync()




class mydict(dict):
	"""
	an extension of the dictionary class. It has the extra function of being able to sort keys according to their
	dictionary value
	"""
	def sorted_keys(self,reverse=False):
		"""
		a function that returns the list of allkeys sorted according to their value.
		Important for histograms (among other things)
		"""
		aux=[(self[k],k)for k in self.keys()]
		aux.sort()
		if reverse:aux.reverse()
		return [k for v, k in aux]		  


##################################################################################################################################################################################
######################## CHEMICAL ORGANIZATION THEORY ############################################################################################################################
##################################################################################################################################################################################


def IsOrganizationSlow(molecules):
	"""
	Returns True if the set of molecules is an organization
	"""
	mol=set(molecules)
	gen=ReactSets(mol,mol)
	if gen==mol:
		return True
	return False







def ExpandOrganizationBy1Molecule(O,molecule):
	if molecule in O:
		return O
	O_New=frozenset([OrgLibrary.check(O|molecule)])
	return O_New

def ExpandOrganizationBy1Organization(O,O2):
	if 02 in O:
		return O
	if O in O2:
		return O2
	O_New=frozenset([OrgLibrary.check(O|02)])
	return O_New


def GetAllKnownOrg():
	return set(NameOrganization.DB.values())

def GetNKnownOrg():
	return len(NameOrganization.DB)

def FindOrgKnownBetween(Omin,Omax):
	#this checks all the orgs and does not use the known information
	return FindOrgKnownAbove(Omin)&FindOrgKnownBelow(Omax)

def FindOrgKnownAbove(Omin):
	AllOrg=GetAllKnownOrg()
	AboveOrg=set([])
	for o in AllOrg:
		if o>Omin:
			AboveOrg|=set([o])
	return AboveOrg

def FindOrgKnownBelow(Omax):
	AllOrg=GetAllKnownOrg()
	BelowOrg=set([])
	for o in AllOrg:
		if o<Omax:
			BelowOrg|=set([o])
	return BelowOrg

def FindOrgByUnionEtIntersection(Orgs):
	"""Given a set of organizations considers all the possible unions and intersections to find all the possible organizations"""
	NewNewOrgs=set([])
	KnownOrgs=copy.deepcopy(Orgs)
	for h in combinations(Orgs,2):
		#checks only if one is not contained in the other
		NewNewOrgs|=frozenset([OrgLibrary.check(h[0]|h[1])])
		#checks only if one is not contained in the other
		NewNewOrgs|=frozenset([OrgLibrary.check(h[0]&h[1])])
	FoundOrgs=NewNewOrgs
	NewOrgs=NewNewOrgs-KnownOrgs
	while NewOrgs:
		NewNewOrgs=set([])
		for h in combinations(NewOrgs,2):
			#checks only if one is not contained in the other
			NewNewOrgs|=frozenset([OrgLibrary.check(h[0]|h[1])])
			#checks only if one is not contained in the other
			NewNewOrgs|=frozenset([OrgLibrary.check(h[0]&h[1])])
		for h in NewOrgs:
			for t in KnownOrgs:
				#checks only if one is not contained in the other
				NewNewOrgs|=frozenset([OrgLibrary.check(h|t)])
				#checks only if one is not contained in the other
				NewNewOrgs|=frozenset([OrgLibrary.check(h&t)])
		KnownOrgs|=NewOrgs
		NewOrgs=NewNewOrgs-KnownOrgs#NewOrgs is what we actually found
	KnownOrgs-=Orgs
	return KnownOrgs

def FindOrgByUnionEtIntersectionKnowingSome(NewOrgsOriginal,KnownOrgsOriginal):
	"""Given a set of organizations considers all the possible unions and intersections to find all the possible organizations"""
	KnownOrgs=copy.deepcopy(KnownOrgsOriginal)
	NewOrgs	 =copy.deepcopy(  NewOrgsOriginal)
	KnownOrgs-=NewOrgs
	while NewOrgs:
		print "Organizations known:",GetNKnownOrg(), "number of new organizations in this research",len(NewOrgs),"Discovered Organization",len(KnownOrgs)
		NewNewOrgs=set([])
		for h in combinations(NewOrgs,2):
			testtime("FindOrgByUnionEtIntersectionKnowingSome h h", OverflowError)
			##			  if profiling:
			##					  if time.time()>timelimit:
			##						  print "FindOrgByUnionEtIntersectionKnowingSome h h"
			##						  raise OverflowError
			NewNewOrgs|=frozenset([OrgLibrary.check(h[0]|h[1])])			#checks only if one is not contained in the other
			NewNewOrgs|=frozenset([OrgLibrary.check(h[0]&h[1])])			#checks only if one is not contained in the other
   #		 print ":",
		print
		print "Organizations known:",GetNKnownOrg(),"recently discovered=",len(NewNewOrgs)
		for h in NewOrgs:
			for t in KnownOrgs:
				if profiling:
					if time.time()>timelimit:
						print "FindOrgByUnionEtIntersectionKnowingSome h t"
						raise OverflowError
				NewNewOrgs|=frozenset([OrgLibrary.check(h|t)])
				NewNewOrgs|=frozenset([OrgLibrary.check(h&t)])
				#			 print ".",
		print
		print "Organizations known:",GetNKnownOrg(),"recently discovered=",len(NewNewOrgs)		  
		KnownOrgs|=NewOrgs
		print "Discovered Organization after intersection",len(KnownOrgs)
		NewOrgs=NewNewOrgs-KnownOrgs#NewOrgs is what we actually found
		print "Discovered Organization after subtractions",len(NewOrgs)
		if profiling:
			if time.time()>timelimit:
				print "FindOrgByUnionEtIntersectionKnowingSome"
				raise OverflowError
	KnownOrgs-=KnownOrgsOriginal
	KnownOrgs-=NewOrgsOriginal
	return KnownOrgs

def FindOrgBetween(Omin, Omax,NMol2Test=2):
	ProduceOrgBetween(Omin, Omax,NMol2Test=NMol2Test)
	return FindOrgKnownBetween(Omin,Omax)

def ProduceOrgBetween(Omin, Omax,NMol2Test=2):
	if len(Omax)<=len(Omin)+1:
		return
	AllOrg=GetAllKnownOrg()
	OrganizationsFound=set([Omin,Omax])
	OtherMolecules=Omax-Omin
	SetsToTest=set([])
	L=len(Omax)
	n=NMol2Test#max number of size of sets to test
	for m in range(1,n+1):
		setsM=[p for p in combinations(OtherMolecules,m)]
		for s in setsM:
			SetsToTest|=frozenset([frozenset(s)|Omin])
			if profiling:
				if time.time()>timelimit:
					print "FindOrgBetween making sets"
					raise OverflowError
	#check out all organizations that are in between
	#check out all the generators for each organization that is in between
		#maybe it is better to just check them all without checking the generators
	for s in SetsToTest:
		OrganizationsFound|=frozenset([OrgLibrary.check(s)])
		if profiling:
			if time.time()>timelimit:
				print "FindOrgBetween"
				raise OverflowError
	FindOrgByUnionEtIntersectionKnowingSome(OrganizationsFound-AllOrg,AllOrg)

def CloseOrganizations():
	"""
	Takes all the organizations known, and adds any organization given from intersection or union that we might have missed
	"""
	FindOrgByUnionEtIntersection(GetAllKnownOrg())

   
def CheckAllSets(molecules):
	"""
	Returns the organizations generated by a set of molecules
	"""
	print "checking all organizations below ", molecules
	Omin=frozenset([])
	Omax=frozenset(OrgLibrary.check(molecules))
	OrganizationsFound=frozenset([Omax,Omin])
	OrganizationsFound|=FindOrgBetween(Omin=Omin,Omax=Omax)
	return OrganizationsFound



IsLatticComplete=False
RecentlyDiscoveredOrgs=set([])
#LatticeLibrary=LatticeLibraryDB()

class MapStatesOrganizations(object):
	"a class that holds the data about the organizations that are present in a series of states"
	states=[0]
	nOrganizations=1
	OrganizationTime=numpy.zeros(shape=(1,1))
	indexstate={}
	def __init__(self, states, nOrganizations=1000):
		#	 def ___new___(self,states,nOrganizations=1000):
		print "called init"
		self.states=states
		self.nOrganizations=nOrganizations
		self.OrganizationTime=numpy.zeros(shape=(self.nOrganizations,len(self.states)))
		indexes=range(len(self.states))
		self.indexstate=dict(zip(indexes, self.states))
		self.stateindex=dict(zip(self.states, indexes))
	def setvalue(self,state,orgname,value):
		try:
			index=self.stateindex[state]
		except KeyError:
			print "aaargh ERROR 234100 You need to expand the map"
		self.OrganizationTime[orgname,index]=value
	def getvalue(self,state,orgname):
		index=self.stateindex[state]
		return self.OrganizationTime[orgname,index]
	def getLastPositiveValue(self,orgname,state):
		"studying the organizations we need to find the last time an organization was active"
		index=self.stateindex[state]
		orghistory=self.OrganizationTime[orgname,0:index+1]#maybe is index, maybe is index-1, to be tested
		i=index #again maybe it should start by index-1
		for t in reversed(orghistory):
			if t:
				return self.indexstate[i]
			i-=1
		return "out of range"#a string is bigger than any number, as such this will not appear when later I will ask for the minimum
	def getOrganizationsActive(self,state):
		index=self.stateindex[state]
		organizations=list(self.OrganizationTime[:,index])
		orgnumbered=zip(organizations,range(self.nOrganizations))
		return [t[1] for t in orgnumbered if t[0]>0]					  
	def getSumLastValues(self,orgname,state):
		"studying the organizations we need to find the last time an organization was active"
		index=self.stateindex[state]
		orghistory=self.OrganizationTime[orgname,0:index+1]#maybe is index, maybe is index-1, to be tested
		return sum(orghistory)	
	def getOrganizationColor(self,state,orgname):
		"""takes a set of organizations, and returns a dictionary that associates to each organization a color
		based on when was the organization active the last time"""
		last=self.getLastPositiveValue(orgname,state)
		if last==0:
			v=0
		else:
			try:
				lastIndex=self.stateindex[last]
				deltaIndex=self.stateindex[state]-lastIndex
				v=1.0/deltaIndex
				v=(float(last)-float(self.indexstate[0]))/(state-self.indexstate[0])
			except ZeroDivisionError:
				v=1
				#				 if float(last)-float(self.indexstate[0]):
				#					 v=1
				#				 else:
				#					 v=1
		lightness=linearapproximation(0.3,.7,v)
		saturation=linearapproximation(0.6,1,v) #how far from the white-gray-black axis mundi 
		hue=linearapproximation(0.6,1,v)
		#print "step",state,"organization ",orgname,"last active on time",last, "color=",hue, saturation,lightness
		return "%f %f %f"%(hue, saturation,lightness)
	def getOrganizationsColors(self,state,organizations):
		"""takes a set of organizations, and returns a dictionary that associates to each organization a color
		based on when was the organization active the last time"""
		lastvalues=[self.getLastPositiveValue(o,state) for o in organizations]
		minvalue=min(lastvalues)
		values=[]
		for lastvalue in lastvalues:
			try:
				vtemp=(float(lastvalue)-float(minvalue))/(state-minvalue)
			except ValueError:
				vtemp="out of range"
			except ZeroDivisionError:
				vtemp=1.0
			values.append(vtemp)
		return [findcolor(v) for v in values]


def findcolor(v):
		try:
			lightness=linearapproximation(0.3,.7,v)
			saturation=linearapproximation(0.6,1,v) #how far from the white-gray-black axis mundi 
			hue=linearapproximation(0.6,1,v)
		except TypeError:
			lightness=1.0
			saturation=0.0
			hue=0.0
		return "%f %f %f"%(hue, saturation,lightness)





###########################################
######## Cell #############################
###########################################

class Cell(object):
	"""
	A Cell class, has the list of molecules inside as 'molecules'. The same list is also present as a mydict dictionary called moleculesdict.
	It has the number of possible reactions as NReactions (which is =-1 if the number has not been calculated),
	and the proportion of reaction that are feasible as NReactionsProp."""
	molecules=[] #a list of all the molecules. len(molecules) will give the size of it. (I chose to use a list instead of a set because it is easier to have multiple elements with the same value
	moleculesdict=mydict()
	NReactions=-1
	NReactionsProp=-0.1
	def __init__(self,molecules=[]):
		"""
		Creates a cell. It takes a list of molecules. If the list is empty it will create one by making a cell of starting_cell_size number of default molecules.
		"""
		self.molecules=[] #a list of all the molecules. len(molecules) will give the size of it. (I chose to use a list instead of a set because it is easier to have multiple elements with the same value
		self.moleculesdict=mydict()
		self.NReactions=-1
		if len(molecules):
			for t in molecules:	 self.addmoleculefromname(t)
		else:
			for t in range(starting_cell_size):
								M=Molecule()
								self.addmolecule(M)
	def CalculateNReactions(self):
		"""
		calculates how many reactions CAN happen in a cell.
		It returns a float given by NReactions divided by number of possible reactions
		"""
		self.NReactions=0
		for n in self.moleculesdict.keys():
			for h in self.moleculesdict.keys():
				if reactionDB.check(n,h) in ReactionsAllowed: #is m produced by n?
					self.NReactions+=self.moleculesdict[n]*self.moleculesdict[h]
		self.NReactionsProp=float(self.NReactions)/(len(self.molecules)**2)
	def derivative(self,m):
		"""
		Given a molecule, calculates the first derivative of the molecule.
		In the context of the cell. It is calculated by seeing how many
		molecules would reproduce a particular molecule.
		"""
		tempP=tempN=0
		for n in self.moleculesdict.keys():
			if reactionDB.check(m,n) in ReactionsAllowed: #is m produced by n?
				tempP+=float(self.moleculesdict[n])/len(self.molecules)
		tempP*=float(self.moleculesdict[m])/len(self.molecules) #and up to here is the positive component.
		if (self.NReactions==-1):
			self.CalculateNReactions()
		tempN=(float(self.moleculesdict[m])/len(self.molecules))*self.NReactionsProp
		return tempP-tempN
	def addmoleculefromname(self,m): 
		""" adds a molecule named m to the cell"""
		self.molecules.append(m)
		try:
			self.moleculesdict[m]+=1
		except KeyError:
			self.moleculesdict[m]=1
		#moleculelibrary.check(m)
	def eliminatemoleculefromname(self,m):
		"""eliminates a molecule named m from the cell"""
		self.molecules.remove(m)
		self.moleculesdict[m]-=1
		if self.moleculesdict[m]==0:
			self.moleculesdict.pop(m)
	def addmolecule(self,M):
		""" adds a molecule M to the cell."""
		self.molecules.append(M.name)
		try:
			self.moleculesdict[M.name]+=1
		except KeyError:
			self.moleculesdict[M.name]=1
		if str(M.name) not in moleculelibrary.DB:
			moleculelibrary.DB[str(M.name)]=M
	def react(self):
		"""picks two molecules, checks if they can react
		finds out the resulting molecule considering the mutations,
		adds it to the Cell,
		and returns true if the reaction happens and false otherwise"""
		global ReactionsAllowed
		reactants=random.sample(self.molecules,3)
		result=react(reactants[0],reactants[1])
		if result:
			mother=moleculelibrary.check(result)
			value,size=mother.reproduce()
			name=valuesize2name(value,size)
			self.addmoleculefromname(name)
			self.eliminatemoleculefromname(reactants[2])			
			return True
		return False
	def biosyntheticphase(self):
		"""picks one molecule, checks if it can absorb the photons
		finds out the resulting molecule considering the mutations,
		adds it to the Cell,
		and returns true if the reaction happens and false otherwise"""
		global ReactionsAllowed
		reactants=random.sample(self.molecules,2)
		probability=photoDB.check(reactants[0])
		if random.random() < probability:
			mother=moleculelibrary[reactants[0]]
			value,size=mother.reproduce()
			name=valuesize2name(value,size)
			self.addmoleculefromname(name)
			self.eliminatemoleculefromname(reactants[1])			
			return True
		return False
	def __repr__(self): return 'Cell('+`self.molecules`+')'
	def __str__(self):
		self.molecules.sort()
		return "['"+str(self.molecules)+"']"
	def prettywrite(self, numbermolecules=12, matrixnumber=1,matrixdecimal=0,PopType=1,showderivative=0 ,matrixplus=0,matrixmutation=0,showorg=0):
		"""
		this function is the key function to represent the state of the cell.
		It returns a string with the state of the cell nicely edited.
		It takes many variables:
			numbermolecules=1...MAX_INT describes how many molecules of the cell should be represented (from the most common onward)
			matrixnumber=0/1/2 describes if the reaction types should be represented 
				0: no,
				1: as a Reaction class Number-0,1,2,...,9-,
				2: as a Reaction CLass symbol.
			matrixdecimal=0/1/2/3/4 describes the way in which the molecules are going to be represented: 
				0: no name at all,
				1: by the name,
				2: by the primary string,
				3: by the secondary string,
				4: by both the primary string and the secondary string.
			PopType=0/1/2 describes the way the population size of each molecule is represented:
				0: none at all
				1: by number
				2: by a series of dots	ex 19=:::::::::.			
			showderivative=0/1	will the first derivative of the function be shown:
				0: no
				1: yes
			matrixplus=0/1/2 independently of the relation between varios molecules, depending on the options.py file some reactions will be made illegal.
				This shows which reactions are possible. Various representations of this info are possible
				0: no representation whatsoever
				1: a '+' sign will appear every time a reaction is possible
				2: given two molecules a and b,
					if ab is reactive and so is ba, then a '*' appears,
					if ab is reactive and ba is not, then a '+' appears at 'ab', and a '-'sign appears at 'ba'
					if neither 'ab' nor 'ba' is reactive nothing appears
			matrixmutation=0/1/2/3 Shows the Lev. mutation distance between two molecules.		  
				0: no representation whatsoever
				1: given two molecules
					at distance 0 show "id" (for identity)
					at distance 1 show "°"
					at distance 2 show "'"
					at distance higher than 2 show nothing
				2: given two molecules
					writes the matrix distance as numbers
				3: given two molecules
					writes the matrix distance as a set of dots, one dot for each mutation necessary to go from one to the other.
			"""
		result=''
		print "ciao",
		keys=self.moleculesdict.sorted_keys(True) #set this to True or False decides the direction of the picture
		if showorg%2==0:
			if showderivative%2:
				generatedmolecules=GeneratedMolecules(keys[:numbermolecules])
			for t in keys[:numbermolecules]:
				if matrixmutation or matrixnumber or matrixplus:
					for h in keys[:numbermolecules]:
						if matrixmutation:						ld=mutationsDB.check(h,t)
						if matrixmutation==0:					result+=""
						elif matrixmutation==1:
							if not ld:						result+="i"
							elif ld==1:						result+="°"
							elif 1<ld<3:					result+="'"
							else:							result+=" "
						elif matrixmutation==2:
							if ld:	result+="% 3s"%ld
							else:	result+=" "
						elif matrixmutation==3:
							if ld:
									resttemp="."+":"*((ld-1)/2)+"."*((ld-1)%2)
									result+=resttemp
							else:	result+=" "
						if	 matrixnumber==2:
							try:
								result+="("+ShowRCN[reactionDB.check(t,h)]	+")"
							except NameError:
								result+=""
						elif matrixnumber==1:
							try:
								result+="("+ShowRC[reactionDB.check(t,h)]  +")"
							except NameError:
								result+=""	  
						elif matrixnumber==0:			result+=""
						result+=writereact(t,h,matrixplus,keys,numbermolecules)
				if matrixplus:
					try:
						result+="%+5.7f\t"%photoDB.check(t)
					except NameError:
						result+=" "	   
						#			 if showderivative==1:
						#				 result+="%+5.7f\t"%self.derivative(t)
				if PopType==1:
					result+="% 5s\t"%self.moleculesdict[t]
				elif PopType==2:
					asterisk=':'*(self.moleculesdict[t]/2)
					asterisk+='.'*(self.moleculesdict[t]%2)
					result+="%s	 "%asterisk
				if showderivative%3:
					try:
						result+="(%i)\t"%generatedmolecules[t]
					except KeyError:
						result+="\t"
				if not matrixdecimal:
					result+="%s"%os.linesep
				elif matrixdecimal==1:
					result+="% 15s"%t+"%s"%os.linesep
				elif matrixdecimal==2:
					result+="%- 50s"%moleculelibrary.check(t).primarystring+"%s"%os.linesep
				elif matrixdecimal==3:
					result+="%- 25s"%moleculelibrary.check(t).secondarystring+"%s"%os.linesep
				elif matrixdecimal==4:
					result+="%- 50s"%moleculelibrary.check(t).primarystring+"	"
					result+="%- 25s"%moleculelibrary.check(t).secondarystring+"%s"%os.linesep
				if numbermolecules<15:
					result+="%s"%os.linesep
		if showorg%2==1:
			Os=LatticeLibrary.check(keys[:numbermolecules])
			#Os=CheckAllSets(keys[:numbermolecules])
			result+="ORGANIZATIONS%s"%(os.linesep)
			Organizations=mydict()
			maxlength=0
			for O in Os:				
				L_O=len(O)
				try:
					Organizations[O]=self.PopulationSet(O)/L_O
				except ZeroDivisionError:
					Organizations[O]=0
				maxlength=max(maxlength,L_O)
			for O in Organizations.sorted_keys(True):
				MoleculesInO=mydict()
				for m in O:
					try:
						MoleculesInO[m]=self.moleculesdict[m]
					except KeyError:
						MoleculesInO[m]=0
				result+="\t"*(1+maxlength-len(O))
				for m in MoleculesInO.sorted_keys(False):
					result+="%s\t"%m
				result+="\t(%s)"%(self.PopulationSet(O))
				result+="\t"*(maxlength/2)
				result+="%s"%(os.linesep)
		return result
	def AddOrganizationsThere(self,step,maporganizations=None):
		"""This function takes the step and the molecules and writes the lattice.
		If it receives only the step it finds the molecules and writes the lattice and draws the diagram,
		if it receives only the molecules it uses the molecules to draw the graph, but does not plot it as a png;
		and if it receives both it draws the graph of the molecules and then stores it as the """
		keys=self.moleculesdict.sorted_keys(True) #set this to True or False decides the direction of the picture
		global Org2FillColor
		Org2Color={}
		Onew=set([])
		#		 Color="%f %f %f"%(random.random(),random.random(),random.random())
		#		 for threashold in range(1,len(keys)+1):
		#			 molecules=keys[:threashold]
		#			 O=OrgLibrary.check(molecules)
		##			  Org2Color[O]="red"
		##			  Org2Color[O]=".1 .1 .99"
		#			 Org2FillColor[O]=Color
		#			 Org2Color[O]="yellow"
		#			 Onew|=set([O])
		#			 molecules=keys[:threashold]
		O=OrgLibrary.check(keys)
		AllOrg=GetAllKnownOrg()
		##		  global RecentlyDiscoveredOrgs
		##		  if RecentlyDiscoveredOrgs:
		##			  Onew=RecentlyDiscoveredOrgs
		##			  AllOrg|=FindOrgByUnionEtIntersectionKnowingSome(Onew,AllOrg)
		##			  RecentlyDiscoveredOrgs=set([])
		##		  AllOrg=GetAllKnownOrg()
		Os=FindOrgKnownBelow(O)|set([O])
		O1=set([O])
		for Organization in Os:
			Org2Color[Organization]="yellow"
		Org2Color[O]="green"
		Org2ColorOnlyO={}
		for Organization in AllOrg:
			Org2ColorOnlyO[Organization]="yellow"
		for Organization in O1:
			Org2ColorOnlyO[Organization]="green"
			#		 Org2ColorOnlyO[O]="green"
		if( maporganizations):
			global OrganizationName
			for Organization in Os:
				oname=OrganizationName.check(Organization)
				maporganizations.setvalue(step,oname,1)
			AllOrgName=[OrganizationName.check(O) for O in AllOrg]
			Colors=maporganizations.getOrganizationsColors(step,AllOrgName)
			Org=[NameOrganization.check(Oname) for Oname in AllOrgName]
			Org2FillColor=dict(zip(Org,Colors))
		WriteLattice(Os,'.',"L_Organizations")
		WriteLattice(AllOrg,'.',"L_AllOrganizations",Org2Color=Org2Color)
		Org2FillColor={}
		WriteLattice(AllOrg,'.',"L_AllOrganizations2",sublattice=O1,Org2Color=Org2ColorOnlyO)
		print "writing the graph"
		#		 print check_call(["dot", "-o./study2/%s.gif"%step, "-Tgif", "./L_AllOrganizations.dot"])
		print check_call(["dot", "-o./study3/%s.jpg"%step, "-Tjpg", "./L_AllOrganizations.dot"])
		print check_call(["dot", "-o./study4/%s.jpg"%step, "-Tjpg", "./L_AllOrganizations2.dot"])
		print "graph drawn"
	def WriteLatticesKnown(self,threashold=0,step=-1,mapStatesOrganizations=None):
		"""This function takes the step and the molecules and writes the lattice.
		If it receives only the step it finds the molecules and writes the lattice and draws the diagram,
		if it receives only the molecules it uses the molecules to draw the graph, but does not plot it as a png;
		and if it receives both it draws the graph of the molecules and then stores it as the """
		keys=self.moleculesdict.sorted_keys(True) #set this to True or False decides the direction of the picture
		if threashold:
			molecules=keys[:threashold]
		else:
			molecules=keys
		if step==-1:
			global stepstudied
			step=stepstudied
		O=OrgLibrary.check(molecules)
		O1=frozenset([O])
		Os=FindOrgKnownBelow(O)|set([O])
		#		  if(mapStatesOrganizations):
		#			 for o in Os:
		#				 print o
		#				 mapStatesOrganizations.setvalue(step,o,1)					  
		AllOrg=GetAllKnownOrg()
		WriteLattice(Os,'.',"L_Organizations")
		print "Os=",Os
		WriteLattice(AllOrg,'.',"L_AllOrganizations",Os)
		print "O1=",O1
		WriteLattice(AllOrg,'.',"L_AllOrganizations2",O1)
		print "writing the graph"
		print check_call(["dot", "-o%s.png"%step, "-Tpng", "./L_AllOrganizations.dot"])
		print "graph drawn"		
	def WriteLattices(self,threashold=0,step=-1):
		"""This function takes the step and the molecules and writes the lattice.
		If it receives only the step it finds the molecules and writes the lattice and draws the diagram,
		if it receives only the molecules it uses the molecules to draw the graph, but does not plot it as a png;
		and if it receives both it draws the graph of the molecules and then stores it as the """
		keys=self.moleculesdict.sorted_keys(True) #set this to True or False decides the direction of the picture
		if threashold:
			molecules=keys[:threashold]
		else:
			molecules=keys
		if step==-1:
			global stepstudied
			step=stepstudied
		Os=LatticeLibrary.check(molecules)		  
		AllOrg=GetAllKnownOrg()
		##		  global IsLatticComplete
		##		  if not IsLatticComplete:
		##			  AllOrg|=FindOrgByUnionEtIntersectionKnowingSome(Os,AllOrg)
		##			  IsLatticComplete=True
		WriteLattice(Os,'.',"L_Organizations")
		WriteLattice(AllOrg,'.',"L_AllOrganizations",Os)
		print "just written L_AltOrg",threashold,step
	def reproduce(self):
		"""a function not really tested yet. It splits the cell into two cells, stores one of them on this cell, and returns the other."""
		metalist=splitlist(self.molecules)
		self.molecules=metalist[0]
		return Cell(metalist[1])
	def calculatebest(self,bests):
		keys=self.moleculesdict.sorted_keys(True) #set this to True or False decides the direction of the picture
		for h in range(1,11):
			try:
				bests[h]|=set(keys[:h])
			except KeyError:
				bests[h]=set(keys[:h])
		for t in range(10,1,-1):
			bests[t]=bests[t]-bests[t-1]
		return bests
	def storemolecules(self,allmolecules,step):
		"""takes all the molecules present in the cell and stores them in a single variable allmolecules"""
		for mol in self.moleculesdict.keys():
			try:
				allmolecules[mol].append((step,self.moleculesdict[mol]))
			except KeyError:
				allmolecules[mol]=[(step,self.moleculesdict[mol])]
	def PopulationSet(self,molecules):
		pop=0
		for m in molecules:
			try:
				pop+=self.moleculesdict[m]
			except KeyError:
				pass
		return pop

def writemolecules(AllMolecules,directory):#not the cell call
		"""
		Takes the integer steps, and the string directory, then writes the number of molecules present. In a series of files.
		One file for each molecule so that eventually executing this file you will have the dictionary of all the times/population size for a certain molecule.
		Very useful to draw the molecule graph.
		The files are written as /molecules/mol_%s.py. where %s represents the molecule name. If the variable 'check' is True,
		and if the value was already present in the file it also checks that the value that we were going to insert was the same that the value that we have found.
		Good to make sure that different runs donw by picking up  the random seed and the state of a system do in fact run in the same way
		"""
		f=open("%s%smolecules.pickle"%(directory,os.sep),'w')
		p = cPickle.Pickler(f)	   # Send pickled data to file f
		p.dump(AllMolecules)
		f.close()


def readmolecules(directory):
		f=open("%s%smolecules.pickle"%(directory,os.sep),'r')
		u = cPickle.Unpickler(f)
		AllMolecules = u.load()
		f.close()
		return AllMolecules		

def splitlist(basiclist, numberoflists=2):
	"""
	splits a list in numberoflists lists and returns a list of lists.
	Each element in the mother list has the same probability to fall in every daughter list
	"""
	newlist=[(t,random.randint(0,numberoflists-1)) for t in basiclist]
	metalist=[[] for t in range(numberoflists)]
	for t in newlist:
		metalist[t[1]].append(t[0])
	return metalist

def getfirstbest(bests,threashold=-1):
	"""
	finds all the molecules that are among the best 'threashold', and stores them into a set.
	The previous version of the set is passed to the function as the first variable, and the union of the two is then returned 
	"""
	if threashold==-1:
		threashold=len(bests)
	b=set([])
	for t in range(1,threashold+1):	 b|=bests[t]
	return b		

def writebest(bests,directory="."):
	"""
	writes the best molecules in a file called bests.py.
	"""
	buf=""
	buf+="bests="+`bests`+"%s"%os.linesep
	f=open("%s%sbests.py"%(directory,os.sep),'w')
	f.write(buf)
	f.close()

def readbest(directory="."):
	"""reads the best molecule files, and make sure that each molecule read is present in the molecule library"""
	f=open("%s%sbests.py"%(directory,os.sep));	 exec f; f.close()
	for t in range(10,1,-1):
		for m in bests[t]:
			M=moleculelibrary.check(m)
	return bests

def getdiversity(step):
	print step
	stepwrittenwell="%10s"%step
	ThisState=AllStates.get(stepwrittenwell, default=None)
	moleculespresent=set(ThisState)
	numbermoleculespresent=len(moleculespresent)
	return numbermoleculespresent
	#return len(set(AllStates["%10s"%step]))
	#	 return len(set(AllStates["%s"%step]))

def getdiversitylist():
	#	val=AllStates.values()
	#	setsmol=[set(v)for v in val]
	#return setsmol
	return 1

def getstates():
	"""gets the list of all the states that has been stores. Sortes them and returns them as a list. Useful to pick the info that is present.
	"""
	b=AllStates.keys()
	b.sort()
	c=[int(d) for d in b]
	return c

def setstate(step):
	print "setting state=",step
	b=Cell(AllStates["%10s"%step])
	random.setstate(AllRndSeeds["%10s"%step])
	return b


#def calculatestate(stepto,stepfrom=0,ReactionsAllowedLocal):
def calculatestate(stepto,stepfrom=0):
	"""
	takes the set of allowed reactions,
	the step we need to claculate [and the step from which we want to start our calculation]
	loads the step we wants to start the claculation from, and then runs the experiement up to the step we are looking for.
	"""
	b=setstate(stepfrom)
	step=stepfrom
	while(step<stepto):#	open("%s/states/state_%s.py"%(directory,step),'w').write("b"+b.writecell())
		step+=1
		b.react()#		  b.react(ReactionsAllowedLocal)
	return b 

def GeneratedMolecules(molecules):
	"""
	Returns the list of all the molecules generated by a set of molecules, with their relative multiplicity
	"""
	results=mydict()
	for m in molecules:
		for n in molecules:
			res=react(m,n)
			if res!= None:
				try:
					results[res]+=1
				except KeyError:
					results[res]=1
	return results





#########################################################################################################
#########################################################################################################
#########################################################################################################

def Org1intersectionOrg2(A,B):
	"""
	Given two organizations (by name), A and B, returns the intersection of the two organizations. If necessary calculating it
	"""
	AiB= IntersectionRelOrg.check(frozenset([A,B]))
	if AiB=="none":
		OrgA=NameOrganization.check(A)
		OrgB=NameOrganization.check(B)
		#		print "OrgA=",OrgA, type(OrgA)
		#		print "OrgB=",OrgB, type(OrgB)		
		OrgAiOrgB = OrgLibrary.check(OrgA & OrgB)
		AiB=OrganizationName.check(OrgAiOrgB)
		IntersectionRelOrg.add(frozenset([A,B]),AiB)
	return AiB

def RelOrg1Org2(A,B):
	"""
	Given two organizations, A and B, returns true if A < B,
	This is equivalent to check if (A intersection B) = A
	"""
	AiB=IntersectionRelOrg.check(frozenset([A,B]))
	if AiB!="none":
		if A==B: return u"="
		if AiB==A: return u"/"
		if AiB==B: return u"\\"
		return u"X"
	else:
		if A==B:
			IntersectionRelOrg.add(frozenset([A,B]),A)
			return u"="
		OrgA=NameOrganization.check(A)
		OrgB=NameOrganization.check(B)
		OrgAiOrgB = OrgA & OrgB
		if OrgAiOrgB==OrgA:
			IntersectionRelOrg.add(frozenset([A,B]),A)
			return u"/"
		if OrgAiOrgB==OrgB:
			IntersectionRelOrg.add(frozenset([A,B]),B)
			return u"\\"
		return u"X"

def FindLostOrganizationByIntersectio(A):
	"""
	sometimes the DB of names gets corrupted. 
	Instead of throwing away all the work we try to recover the data 
	by looking for a couple of organization whose intersection form that organization. 
	Then we recalculate the intersection and use that
	"""
	

def IsOrg1inOrg2(A,B):
	"""
	Given two organizations, A and B, returns true if A < B,
	This is equivalent to check if (A intersection B) = A
	"""
	AiB=IntersectionRelOrg.check(frozenset([A,B]))
	if AiB!="none":
		return AiB==A
	else:
		OrgA=NameOrganization.check(A)
		OrgB=NameOrganization.check(B)
		OrgAiOrgB = OrgA & OrgB
		if OrgAiOrgB==OrgA:
			IntersectionRelOrg.add(frozenset([A,B]),A)
			return 1
		return 0

def MinOrg(SetOrg):
	"""
	given a set of organizations, passed by name,
	returns the list of all organizations that don't have another organization in the list
	smaller than them.
	essentially is a pareto front with the dominance being that A dominates B is A<B
	"""
	Result=SetOrg
	for l in SetOrg:
		for h in ListOrg:
			if h==l: continue
			if IsOrg1inOrg2(h,l):
				Result.remove(l)
				break
	return Result

def MaxOrg(SetOrg):
	"""
	given a set of organizations, passed by name,
	returns the list of all organizations that don't have another organization in the list
	BIGGER than them.
	essentially is a pareto front with the dominance being that A dominates B is A>B
	"""
	Result=SetOrg
	for l in SetOrg:
		for h in ListOrg:
			if h==l: continue
			if IsOrg1inOrg2(l,h):
				Result.remove(l)
				break
	return Result





#########################################################################################################
###########################################Calc Organizations##############################################################
#########################################################################################################



#TableIntersection={} #Given two organizations (A,B), stored by name, returns the union of them
#TableUnion={} #Given two organizations (A,B), stored by name, returns the union of them

#S

#-TableRelations={} #Given two organizations (A,B), returns the relation between them:
##3 =  A greater than B		-3 = B greater than A
##2 = A possibly covers B		-2 = B possibly covers A
##1 = A surely covers B		-1 = B surely covers A
##“I” = Identity (the diagonal)		“<>” = A and B are uncomparable
##			0 = no information



UnsolvedIntersections=set([])
UnsolvedUnions=set([])






#Generators={} #the minimum sets that generate something

timeresult=[]
timestart=time.time()
starttime=0.0




def UpdateTime(HappeningNow=()):
	global timeresult
	global timestart
	global OrganizationsStudied
	informations=(time.time()-timestart,len(OrganizationsStudied))+tuple(HappeningNow)
	timeresult.append(informations)

def PrintTime():
	global timeresult
	print
	print
	for t in timeresult:	
		for r in t:
			print r,"\t",
		print
	print
	



def FindAllTriangles(Orgs=None):
	"""COST:|Orgs|*|OCovered(Org)|*(|SETUNION(OMaybeCovered+OCovered)|+|SETINTERSECTION(OMaybeCovered*2)|)
	"""
	if Orgs==None: 	Orgs=OrganizationsStudied
	for Org in Orgs:
		FindTriangles(Org)

def FindTriangles(Org):
	"""COST:|OCovered(Org)|*(|SETUNION(OMaybeCovered+OCovered)|+|SETINTERSECTION(OMaybeCovered*2)|)"""
	for t in OCovered[Org]:
		Triangles=(OMaybeCovered[t]|OCovered[t])&OMaybeCovered[Org] 
		for s in Triangles:
			WriteLatticeCalculated(set([Org,t,s]),name="Triangle%s_%s_%s"%(Org,t,s))
			print "Triangle Found: %s >- %s >- %s"%(Org,t,s)
			#	for t in OrganizationsCovered:
			#		(OrganizationsMaybeCoveredT,OrganizationsCoveredT,OrganizationsCoveringT,OrganizationsMaybeCoveringT)=OrganizationsNear(t)
			#		Triangles=(OrganizationsMaybeCoveredT|OrganizationsCoveredT)&OrganizationsMaybeCovered
			#		for s in Triangles:
			#			WriteLatticeCalculated(set([Org,t,s]),name="Triangle%s_%s_%s"%(Org,t,s))
			#			print "Triangle Found: %s >- %s >- %s"%(Org,t,s)
			#	print "Triangles on %s checked"%Org

def WriteLatticeCalculated(S,name="LatticeFoundDefaultName",directory=".",DrawCovered=True,DrawMaybeCovered=True):
	#-global TableRelations
	global NameSize
	buf="digraph G {%s%s"%(os.linesep,os.linesep)
	layers=[]
	r0=""
	for s in S:
		try:
			if NameSize[s]<1:	continue
		except KeyError:		
			NameSize[s]=len(NameOrganization.check(s))				
		r0='{rank=same; '
		r0+="%s "%s
		r0+='"_%s"'%NameSize[s]
		layers.append(NameSize[s])
		
		r0+='};%s'%os.linesep
		buf+=r0
	layers=set(layers)
	layers=list(layers)
	layers.sort()
	layers.reverse()
	for lay in range(0,len(layers)-1):		r0+='"_%s" -> "_%s" [color=white];%s'%(layers[lay],layers[lay+1],os.linesep)
	for lay in range(0,len(layers)):		r0+='"_%s" [color=white];%s'%(layers[lay],os.linesep)
	buf+=r0
	for r in S:
		if DrawCovered:
			for s in OCovered[r]&S:
				r0=''				
				r0+="%s -> "%r
				r0+="%s "%s
				r0+="[color=red]"
				r0+="%s"%(os.linesep)
				buf+=r0
		if DrawMaybeCovered:
			for s in OMaybeCovered[r]&S:
				r0=''				
				r0+="%s -> "%r
				r0+="%s "%s
				r0+="[color=black]"
				r0+="%s"%(os.linesep)
				buf+=r0		
	buf+="}%s%s"%(os.linesep,os.linesep)
	f=open("%s%s%s.dot"%(directory,os.sep,name),'w')
	f.write(buf)
	f.close()

def OrganizationsDifferences(Org1,Org2):
	OrgDetails1=NameOrganization.check(Org1)
	OrgDetails2=NameOrganization.check(Org2)
	return (OrgDetails1-OrgDetails2,OrgDetails2-OrgDetails1)



#def CheckUnions(Orgs1=None,Orgs2=None):
	#	"DO NOT USE THIS"
	#	if Orgs1==None: 	Orgs1=OrganizationsStudied
	#	if Orgs2==None: 	Orgs2=OrganizationsStudied
	#	for OrgA in Orgs1:
	#		for OrgB in Orgs2:
	#			if (OrgA,OrgB) in TableUnion:			
	#				OrgU=Org1UnionOrg2(OrgA,OrgB)
	#				if TableUnion[(OrgA,OrgB)]!=OrgU or TableUnion[(OrgB,OrgA)]!=OrgU: 
	#					print "Org %s U %s = %s but the Table returns %s (with opposite %s)"%(OrgA,OrgB,OrgU,TableUnion[(OrgA,OrgB)],TableUnion[(OrgB,OrgA)])

def CheckSingleRelation(OrgDetails1,OrgDetails2):
	if OrgDetails1>OrgDetails2:		return    1
	elif OrgDetails1<OrgDetails2:	return   -1
	elif OrgDetails1==OrgDetails2:	return  "I"
	else: 							return "<>"

def CheckSingleRelationName(Org1,Org2):
	return CheckSingleRelation(NameOrganization.check(Org1),NameOrganization.check(Org2))

def CheckRelation(Org1,Org2):
	#-TableResult=TableRelations[(Org1,Org2)]
	TableResult=CalculateRelation(Org1,Org2)
	ActualRelation=CheckSingleRelation(NameOrganization.check(Org1),NameOrganization.check(Org2))
	if TableResult in (1,2,3) and ActualRelation!=1:
		print Org1, Org2, TableResult, ActualRelation
		assert 0
	elif TableResult in (-1,-2,-3) and ActualRelation!=-1:
		print Org1, Org2, TableResult, ActualRelation
		assert 0
	elif  TableResult=="I" and ActualRelation!="I":
		print Org1, Org2, TableResult, ActualRelation
		assert 0
	elif  TableResult=="<>" and ActualRelation!="<>":
		print Org1, Org2, TableResult, ActualRelation
		assert 0
		

def CheckSomeRelations(Orgs1,Orgs2):
	for O1 in Orgs1:
		for O2 in Orgs2:
			CheckRelation(O1,O2)

def TestAllRelations(Orgs=OrganizationsStudied):
	result=True
	Buf=""
	for O in Orgs:
		Buf+=CheckRelations(O)
	return Buf

def CheckRelations(Org):
	Buf=""
	OrgDetails=NameOrganization.check(Org)
	for s in OrganizationsStudied:
		sDetails=NameOrganization.check(s)
		rel=CheckSingleRelation(OrgDetails,sDetails)
		relTable=CalculateRelation(Org,s)
		tableRel=CalculateRelation(s,Org)
		
		#-if rel==1 and TableRelations[(Org,s)] not in (1,2,3):
		if rel==1 and relTable not in (1,2,3):
			tableRel=CalculateRelation(s,Org)
			Buf+="Org %s contains %s but the Table returns %s (with opposite %s)%s"%(Org,s,relTable,tableRel,os.linesep)
		#elif rel==-1 and TableRelations[(Org,s)] not in (-1,-2,-3):
		elif rel==-1 and relTable not in (-1,-2,-3):
			Buf+="Org %s is contained in %s but the Table returns %s (with opposite %s)%s"%(Org,s,relTable,tableRel,os.linesep)
		#elif rel=="I" and TableRelations[(Org,s)]!="I":
		elif rel=="I" and relTable!="I":
			Buf+="Org %s is equal to itself (%s) but the Table returns %s (with opposite %s)%s"%(Org,s,relTable,tableRel,os.linesep)
		#elif rel=="<>" and TableRelations[(Org,s)]!="<>":
		elif rel=="<>" and relTable!="<>":
			Buf+="Org %s is incomparable with %s but the Table returns %s (with opposite %s)%s"%(Org,s,relTable,tableRel,os.linesep)
	return Buf
			

def CheckAllRelations(Orgs=OrganizationsStudied):
	Buf=""
	for s in Orgs:	
		Buf+=CheckRelations(s)
	return Buf

def DrawTemporarily(setNames=OrganizationsStudied,filename="orgStudied",directory=".",DrawCovered=True,DrawMaybeCovered=True):
	WriteLatticeCalculated(setNames,name=filename+"_R",directory=directory,DrawCovered=DrawCovered,DrawMaybeCovered=DrawMaybeCovered)

def CheckStateOrganizations():
	for t in OrganizationsStudied:
		print t, NameOrganization.check(t)



def PrintReactionScheme(S,WriteAll=False):
	"returns all the reactions in a set as a string"
	buf="["
	for m in S:
		for n in S:
			mn=react(m,n)
			if mn or WriteAll:
				buf+="%s+%s->%s;"%(m,n,mn)
	buf+="]"		
	return buf

def AllMolecules(SetOfSets):
	"Given a set of sets of molecules returns all the molecules contained"
	result=set([])
	for t in SetOfSets:
		result|=t
	return result

def AllMoleculesInOrgs(Orgs):
	Mols=set([])
	for o in Orgs:
		Mols|=NameOrganization.check(o)
	return Mols
		

def AreAlive(SS):
	"Returns False if all the molecules inside have no reactions among them, if there is at least a single reaction returns true"
	for S in SS:
		if IsAlive(S):
			return True
	return False

def IsAlive(S):
	"Returns False if all the molecules inside have no reactions among them, if there is at least a single reaction returns true"
	for m in S:
		for n in S:
			if react(m,n)!=None:
				return True 
	return False





	






def GetSetsRelationshipAll(Orgs,R,Context=OrganizationsStudied):
	Result=Context
	for O in Orgs:
		Result&=R[O]
	return Result







def CalculateRelation(A,B):
	if 		B in OCovering[A]:			return -1
	elif	B in OMaybeCovering[A]:		return -2
	elif 	B in OAboveNotCovering[A]:	return -3
	elif 	B in OCovered[A]:			return  1
	elif 	B in OMaybeCovered[A]:		return  2
	elif 	B in OBelowNotCovered[A]:	return  3
	elif 	B in OSide[A]:				return "<>"
	elif	A==B:	return "I"
	assert 0





def CheckOrganizationsNextLayer(LayerBefore,NewOrganizations,TopOrg,BottomOrg,BottomOrgDetails,MaxNumberOrg=0):
	"""Given a layer of sets of size n that fail to generate any organization, this function calls the function to generate the next layer to test,
	and then tests each element in that layer. Returns the elements of the next layer that fail to generate something, 
	as well as the new organizations generated, united with the previous organizations generated"""
	NextLayer=set([])
	SetsToTest=NextSetsToTest(LayerBefore)
	print "SetsToTest=",len(SetsToTest),
	if len(SetsToTest):
		print " example: ",list(SetsToTest)[0]
	else: print
	for FrozenSetToTest in SetsToTest:
		SetToTest=FrozenSetToTest
		for o in OMaybeCovered[TopOrg]|OCovered[TopOrg]:
			if SetToTest in NameOrganization.check(o):
				print "set to test ",SetToTest,"is already present in ",o
				print "o=",NameOrganization.check(o)
				continue
		S=set(FrozenSetToTest)|BottomOrgDetails
		OrgS=Organization(S)
		if OrgS==BottomOrgDetails:
			NextLayer|=set([FrozenSetToTest])
			continue
		else:
			NameS=OrganizationName.check(OrgS)
			if NameS not in NewOrganizations:
				NewOrganizations|=set([NameS])
				print DescribeSituation(),
				ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections=AddOrganizationOfficially(NameS,TopOrg,BottomOrg)
				if ToBeSolvedByHandUnions:					ToBeSolvedByHandUnionsDict[NameS]=ToBeSolvedByHandUnions
				if ToBeSolvedByHandIntersections:			ToBeSolvedByHandIntersectionsDict[NameS]=ToBeSolvedByHandIntersections				
				print DescribeSituation(),
				sys.stdout.flush()
				while ToBeSolvedByHandUnionsDict or ToBeSolvedByHandIntersectionsDict:
					UpdateTime()
					PrintTime()
					
					keyslength=[(k,NameSize[k]) for k in ToBeSolvedByHandUnionsDict.keys()]
					keyslength.sort(sortpairs, reverse=False)
					for k,s in keyslength:
						NewOrganizations|=SolveByHandUnion(k)
						if MaxNumberOrg and len(NewOrganizations)>MaxNumberOrg:
							UpdateTime()
							PrintTime()
								
							return	set([]),NewOrganizations
					keyslength=[(k,NameSize[k]) for k in ToBeSolvedByHandIntersectionsDict.keys()]
					keyslength.sort(sortpairs, reverse=True)
					for k,s in keyslength:
						NewOrganizations|=SolveByHandIntersection(k)
						if MaxNumberOrg and len(NewOrganizations)>MaxNumberOrg:	
							UpdateTime()
							PrintTime()
							
							return	set([]),NewOrganizations
	return NextLayer,NewOrganizations

def NextSetsToTest(PreviousGoods):
	"""Given a set of elements, none of those that makes an organization,
	their combination might make an organization. 
	The combination in pairs of single elements is given by:set([frozenset(t) for t in itertools.combinations(PreviousGoods, 2)])
	If instead I am given sets of size n, I need to build sets of size n+1, such that each subset of size n is among the sets I am testing.
	Because if a combination is not contained in the set it means that it has already been tested. So each combination that contains that 
	would either reproduce the same or produce something bigger. In which case we don't need to test it here.
	The function takes either a set of elements, or a set of frozensets of elements all of the same size n. And returns a set of sets of elements.
	Each set contained having of size n+1.
	"""
	if len(PreviousGoods)==0: 					return set([])
	try:
		LenAim=len(list(PreviousGoods)[0])+1
	except TypeError: #LenAim=2; But PreviousGoods are here int and not sets
		return set([frozenset(t) for t in itertools.combinations(PreviousGoods, 2)])
	result=set([])
	LenAimIntersection=LenAim-3
	for i in PreviousGoods:
		for j in PreviousGoods:
			ToTest=i|j
			if len(ToTest)!=LenAim:		continue			#next j
			SymmetricDiff=i^j
			IntersectionIJ=i&j			
			for h in itertools.combinations(IntersectionIJ, LenAimIntersection):
				if frozenset(h)|SymmetricDiff in PreviousGoods:							#					print "testing:",ToTest," just found ",frozenset(h)|SymmetricDiff," inside previous layer"	
					continue
				else:																	#					print "testing:",ToTest," Could not find ",frozenset(h)|SymmetricDiff," in previous layer"																
					break 				#next j
			else: result|=set([ToTest])
	return result






def FindNewOrganizationInBetween(TopOrg,BottomOrg,MaxNumberOrg=None):
	"""Given two organizations TopOrg, BottomOrg, with BottomOrg contained in TopOrg,
	finds all the organizations generated between toporg and bottomorg... well eventually
	"""
	TopOrgDetails=NameOrganization.check(TopOrg)
	BottomOrgDetails=NameOrganization.check(BottomOrg)
	ToBeCheckedMolecules=TopOrgDetails-BottomOrgDetails
	MoleculesTodo=len(ToBeCheckedMolecules) #actually they are molecules, not atoms :-)
	NewOrganizations=set([])
	#ThisLayerOrganizations=set([])
	print "Looking for the organisations above ",BottomOrg," and below ",TopOrg,"; Molecules to be checked: ",MoleculesTodo,": ",ToBeCheckedMolecules
	Implosive=[]	#those are the downward molecules
	OriginalOrganizationsFound=0
	UpdateTime()
	PrintTime()
	MolsCovered=set([])
	for a in ToBeCheckedMolecules:
		#DrawTemporarily(filename="orgStudied%s"%MoleculesTodo,DrawMaybeCovered=False)
		MoleculesTodo-=1
		ThisLayerOrganizations=set([])
		if a in MolsCovered:
			print "Molecule",a, "was already in Molecules Covered (MolsCovered)"
			continue
		print "Looking for the organisations above ",BottomOrg," and below ",TopOrg,"; Molecule to be checked: ",a
		S=set([a])|BottomOrgDetails
		OrgS=Organization(S)
		if OrgS==BottomOrgDetails:
			Implosive.append(a)
			continue
		else:
			NameS=OrganizationName.check(OrgS)						#			ExplosiveSets|=set([frozenset([a])])
			if NameS not in NewOrganizations:
				OriginalOrganizationsFound+=1
				NewOrganizations|=set([NameS])
				#ThisLayerOrganizations|=set([NameS])
				ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections=AddOrganizationOfficially(NameS,TopOrg,BottomOrg)
				if ToBeSolvedByHandUnions:					ToBeSolvedByHandUnionsDict[NameS]=ToBeSolvedByHandUnions
				if ToBeSolvedByHandIntersections:			ToBeSolvedByHandIntersectionsDict[NameS]=ToBeSolvedByHandIntersections
				print DescribeSituation(),
				LogThis()
				sys.stdout.flush()
				
				#we found an organisation, so now we complete the lattice
				while ToBeSolvedByHandUnionsDict or ToBeSolvedByHandIntersectionsDict:
					#UpdateTime()
					#PrintTime()
					
					keyslength=[(k,NameSize[k]) for k in ToBeSolvedByHandUnionsDict.keys()]
					keyslength.sort(sortpairs, reverse=False)
					for k,s in keyslength:
						NewOrganizations|=SolveByHandUnion(k)
						if MaxNumberOrg and len(NewOrganizations)>MaxNumberOrg:	
							#UpdateTime()
							#PrintTime()
							
							MolsCovered=AllMoleculesInOrgs(OMaybeCovered[TopOrg]|OCovered[TopOrg])
							print "MolsCovered:",MolsCovered
							#print "Atoms left:",MoleculesTodo
							print "Molecules left:",MoleculesTodo,
							print "actual Molecules to check",len(ToBeCheckedMolecules-MolsCovered)
							
							return	
						LogThis()
					keyslength=[(k,NameSize[k]) for k in ToBeSolvedByHandIntersectionsDict.keys()]
					keyslength.sort(sortpairs, reverse=True)
					for k,s in keyslength:
						NewOrganizations|=SolveByHandIntersection(k)
						if MaxNumberOrg and len(NewOrganizations)>MaxNumberOrg:	
							#UpdateTime()
							#PrintTime()
							MolsCovered=AllMoleculesInOrgs(OMaybeCovered[TopOrg]|OCovered[TopOrg])
							print "MolsCovered:",MolsCovered
							#print "Atoms left:",MoleculesTodo
							print "Molecules left:",MoleculesTodo
							print "actual Molecules to check",len(ToBeCheckedMolecules-MolsCovered)
							
							return	
						LogThis()
					
					LogThis()				
				MolsCovered=AllMoleculesInOrgs(OMaybeCovered[TopOrg]|OCovered[TopOrg])
				print "MolsCovered:",MolsCovered
				print "Molecules left:",MoleculesTodo
				print "actual Molecules to check",len(ToBeCheckedMolecules-MolsCovered)
				#	CheckAllRelations()
				#	for s in NewOrganizations:      #		print "re working on %s",s
				#		AddOrganizationUnionIntersection2(s)
				#	CheckAllRelations()
	MolsCovered=AllMoleculesInOrgs(OMaybeCovered[TopOrg]|OCovered[TopOrg])
	print "MolsCovered:",MolsCovered
	print "Molecules left:",MoleculesTodo
	print "actual molecules to check",len(ToBeCheckedMolecules-MolsCovered)

	ImplosiveSets=set(Implosive)	
	#
	print DescribeSituation()
	#	print "UnsolvedUnions ",len(UnsolvedUnions)," UnsolvedIntersections ",len(UnsolvedIntersections)
	print "I am about to go to the level 2"
	sys.stdout.flush()
	ImplosiveSets,NewOrganizations=CheckOrganizationsNextLayer(ImplosiveSets,NewOrganizations,TopOrg,BottomOrg,BottomOrgDetails,MaxNumberOrg)
	#
	print DescribeSituation()
	print "ImplosiveSets=",len(ImplosiveSets),"Example: ",list(ImplosiveSets)[0],"Is any of those sets Alive?",AreAlive(ImplosiveSets)
	while(AreAlive(ImplosiveSets)):
		ImplosiveSets,NewOrganizations=CheckOrganizationsNextLayer(ImplosiveSets,NewOrganizations,TopOrg,BottomOrg,BottomOrgDetails,MaxNumberOrg)
		print DescribeSituation()
		print "ImplosiveSets=",len(ImplosiveSets),"Example: ",list(ImplosiveSets)[0],"Is any of those sets Alive?",AreAlive(ImplosiveSets)
		sys.stdout.flush()
	#CheckAllRelations()
	#FindAllTriangles()
	#We need to store that this relation has been exhausted. so the relations hould pass from =2 to equal 3 (if any organization has been found)
	#	
	sys.stdout.flush()
	return NewOrganizations
















#from itertools import chain, combinations



	
def CheckAllSubsetsBySize(Molecules):
	"Checks one by one all the subsets of molecules, to find which one are organisations. But then returns them without storing them"
	organisationsFound=set([])
	StartTime=time.clock()
	NOrgFound=0
	for l in range(len(Molecules)):
		SetsToTest=combinations(Molecules, l)
		print "testing size",l
		for s in SetsToTest:
			if IsOrganization(s):
					organisationsFound.add(frozenset(s))
					NowTime=time.clock()
					TimeNeeded=NowTime-StartTime
					NOrgFound+=1
					AverageTime=TimeNeeded/NOrgFound
					print TimeNeeded, " Organisation found: ", NOrgFound, ". Average Time:", AverageTime, " Content:", s
	return organisationsFound




	
def CompareOrgFoundTwoMethods(Orgs1,Orgs2):
	"given two sets of organisations, checks what organisations are in common and what are on one or on the other. It has never been tested so far yet"
	OrgsIntersection=set({})
	Orgs1Minus2=set([])
	Orgs2Minus1=set([])

	for O1 in Orgs1Changeable:
		for O2 in Orgs2:
			if O1==O2:
				OrgsIntersection.add(O1)
				break
	Orgs1Minus2=O1-OrgsIntersection
	Orgs2Minus1=O2-OrgsIntersection
	return Orgs1Minus2, Orgs2Minus1, OrgsIntersection

		

	
#Given two organizations A and B with A>B, And given a set of organizations S1, S2, ..., Sn, I only need to check the elements
#that are not contained in any of the Si (as single elements). Then I can check the pair of elements (p1, p2) with p1 in Si, p2 in Sj, i!=j, and Si intersection Sj = A and Si U Sj=B
#Then I can build up from those. *****HERE I NEED TO CONTINUE*******


def MapLattice(NameTopOrg, NameBottomOrg,MaxNumberOrg=None):
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
		
	#DrawTemporarily()
		
	#-TableRelations[(NameTopOrg, NameBottomOrg)]=2
	#-TableRelations[(NameBottomOrg,NameTopOrg)]=-2
	AddSetRelation(NameTopOrg,NameBottomOrg,(OMaybeCovered,OMaybeCovering))

	#-TableRelations[(NameTopOrg,NameTopOrg)]="I"
	#-TableRelations[(NameBottomOrg,NameBottomOrg)]="I"

	TableUnion[(NameTopOrg,NameBottomOrg)]=NameTopOrg
	#	TableUnion[(NameBottomOrg,NameTopOrg)]=NameTopOrg
	TableUnion[(NameTopOrg,NameTopOrg)]=NameTopOrg
	TableUnion[(NameBottomOrg,NameBottomOrg)]=NameBottomOrg

	#	STableUnion.setvalue((NameTopOrg,   NameBottomOrg),    NameTopOrg)
	#	STableUnion.setvalue((NameTopOrg,      NameTopOrg),    NameTopOrg)
	#	STableUnion.setvalue((NameBottomOrg,NameBottomOrg), NameBottomOrg)
	

	TableIntersection[(NameTopOrg,NameBottomOrg)]=NameBottomOrg
	#	TableIntersection[(NameBottomOrg,NameTopOrg)]=NameBottomOrg
	TableIntersection[(NameBottomOrg,NameBottomOrg)]=NameBottomOrg
	TableIntersection[(NameTopOrg,NameTopOrg)]=NameTopOrg

	OrgInBetween=FindNewOrganizationInBetween(NameTopOrg,NameBottomOrg,MaxNumberOrg)
	print OrgInBetween





	
	#OrgInBetween=FindNewOrganizationsAbove(NameBottomOrg,molecules,NameTopOrg)
	
	
	
	
	#print OrgInBetween


def FindNewOrganizationsAboveByX(BottomOrg,molecules,TopOrg,OrganizationsKnown,SizeSets=2,MaxNumberOrg=None):
	"""made redundant by FindNewOrganizationsAboveByY
	Given an organizations BottomOrg, 
	finds all the organizations generated by bottomorg adding SizeSets molecules at a time.
	and each time it finds an organisation it expands the sub-lattice
	"""
	
	StartTimeX=time.clock()
	
	
	
	BottomOrgDetails=NameOrganization.check(BottomOrg)
	try:
		ToBeCheckedMolecules=DownwardMolecules[BottomOrg][1]
	except KeyError:
		return set([])
	MoleculesTodo=len(ToBeCheckedMolecules)
	NewOrganizations=set([])
	#ThisLayerOrganizations=set([])
	print "Looking for the organisations above ",BottomOrg,"; Molecules to be checked: ",MoleculesTodo,": ",ToBeCheckedMolecules
	print "  UpwardMolecules:",len(  UpwardMolecules[BottomOrg][1]),":",  UpwardMolecules[BottomOrg][1]
	print "SidewardMolecules:",len(SidewardMolecules[BottomOrg][1]),":",SidewardMolecules[BottomOrg][1]
	print "DownwardMolecules:",len(DownwardMolecules[BottomOrg][1]),":",DownwardMolecules[BottomOrg][1]

	if SizeSets not in DownwardMolecules[BottomOrg].keys():         DownwardMolecules[BottomOrg][SizeSets]=set([])
	if SizeSets not in UpwardMolecules[BottomOrg].keys():             UpwardMolecules[BottomOrg][SizeSets]=set([])
	if SizeSets not in SidewardMolecules[BottomOrg].keys():         SidewardMolecules[BottomOrg][SizeSets]=set([])

	print "  UpwardMolecules:[",SizeSets,"]",len(  UpwardMolecules[BottomOrg][SizeSets]),":",  UpwardMolecules[BottomOrg][SizeSets]
	print "SidewardMolecules:[",SizeSets,"]",len(SidewardMolecules[BottomOrg][SizeSets]),":",SidewardMolecules[BottomOrg][SizeSets]
	print "DownwardMolecules:[",SizeSets,"]",len(DownwardMolecules[BottomOrg][SizeSets]),":",DownwardMolecules[BottomOrg][SizeSets]


	#	DeadEnd=set()
	#	DeadEnd=SidewardMolecules[BottomOrg][SizeSets-1]|UpwardMolecules[BottomOrg][SizeSets-1]

	OriginalOrganizationsFound=0

	print "ToBeCheckedMolecules=",ToBeCheckedMolecules
	SetsToTest=findsubsets(ToBeCheckedMolecules,SizeSets)
	
	
	print " ",len(SetsToTest)," Sets to Test: ",SetsToTest


	#	MolsCovered=DownwardMolecules[BottomOrg]|UpwardMolecules[BottomOrg]|SidewardMolecules[BottomOrg]
	#	ToBeCheckedMolecules=ToBeCheckedMolecules-MolsCovered
	for a in SetsToTest:
		MoleculesTodo-=1
		ThisLayerOrganizations=set([])	
		print a
		SetA=set(a)	
		print "Looking for the organisations above ",BottomOrg,"; Set to be checked: ",SetA,
		print "SetA=",SetA
		S=SetA|BottomOrgDetails
		print "S=",S
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


def FindNewOrganizationsAboveBy2(BottomOrg,molecules,TopOrg,OrganizationsKnown,MaxNumberOrg=None):
	"""made redundant by FindNewOrganizationsAboveByY
	Given an organizations BottomOrg, 
	finds all the organizations generated by bottomorg adding two molecule at a time.
	and each time it finds an organisation it expands the sub-lattice
	"""
	
	SizeSets=2
	
	StartTimeX=time.clock()
	
	
	BottomOrgDetails=NameOrganization.check(BottomOrg)
	try:
		ToBeCheckedMolecules=DownwardMolecules[BottomOrg][1]
	except KeyError:
		return set([])
	MoleculesTodo=len(ToBeCheckedMolecules)
	NewOrganizations=set([])
	#ThisLayerOrganizations=set([])
	print "Looking for the organisations above ",BottomOrg,"; Molecules to be checked: ",MoleculesTodo,": ",ToBeCheckedMolecules
	print "  UpwardMolecules:",len(  UpwardMolecules[BottomOrg][1]),":",  UpwardMolecules[BottomOrg][1]
	print "SidewardMolecules:",len(SidewardMolecules[BottomOrg][1]),":",SidewardMolecules[BottomOrg][1]
	print "DownwardMolecules:",len(DownwardMolecules[BottomOrg][1]),":",DownwardMolecules[BottomOrg][1]

	if SizeSets not in DownwardMolecules[BottomOrg].keys():         DownwardMolecules[BottomOrg][SizeSets]=set([])
	if SizeSets not in UpwardMolecules[BottomOrg].keys():             UpwardMolecules[BottomOrg][SizeSets]=set([])
	if SizeSets not in SidewardMolecules[BottomOrg].keys():         SidewardMolecules[BottomOrg][SizeSets]=set([])

	print "  UpwardMolecules:[",SizeSets,"]",len(  UpwardMolecules[BottomOrg][SizeSets]),":",  UpwardMolecules[BottomOrg][SizeSets]
	print "SidewardMolecules:[",SizeSets,"]",len(SidewardMolecules[BottomOrg][SizeSets]),":",SidewardMolecules[BottomOrg][SizeSets]
	print "DownwardMolecules:[",SizeSets,"]",len(DownwardMolecules[BottomOrg][SizeSets]),":",DownwardMolecules[BottomOrg][SizeSets]


	#	DeadEnd=set()
	#	DeadEnd=SidewardMolecules[BottomOrg][SizeSets-1]|UpwardMolecules[BottomOrg][SizeSets-1]

	if MoleculesTodo<SizeSets:  #in the most general case we need to check that there are at least n sets of size n-1 as a requirement for there be a set of size n
		return NewOrganizations
	OriginalOrganizationsFound=0

	print "ToBeCheckedMolecules=",ToBeCheckedMolecules
	SetsToTest=findsubsets(ToBeCheckedMolecules,SizeSets)
	
	
	print " ",len(SetsToTest)," Sets to Test: ",SetsToTest


	#	MolsCovered=DownwardMolecules[BottomOrg]|UpwardMolecules[BottomOrg]|SidewardMolecules[BottomOrg]
	#	ToBeCheckedMolecules=ToBeCheckedMolecules-MolsCovered
	for a in SetsToTest:
		MoleculesTodo-=1
		ThisLayerOrganizations=set([])	
		print a
		SetA=set(a)	
		print "Looking for the organisations above ",BottomOrg,"; Set to be checked: ",SetA,
		print "SetA=",SetA
		S=SetA|BottomOrgDetails
		print "S=",S
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



########################NOT NECESSARY##############



def GetRAll(Orgs,R,Context=OrganizationsStudied):
	"""Returns the Organizations in Context that have a relations R with all of the organizations in Orgs"""
	assert 0
	Result=[]
	for P in Context:
		for O in Orgs:
			if TableRelations[(O,P)] not in R:
				break
		else:
			Result.append(P)
	return set(Result)

def GetRAny(Orgs,R,Context=None):
	"""Returns the Organizations in Context that have a relations R with any of the organizations in Orgs, Essentially it is a Union"""
	assert 0
	if Context==None:	Context=OrganizationsStudied
	Result=[]
	for P in Context:
		for O in Orgs:
			if TableRelations[(O,P)] in R:
				Result.append(P)
	return set(Result)

def AddRelation(Org, Orgs, Relation):
	assert 0
	for O in Orgs:
		TableRelations[(Org,O)]=Relation[0]
		TableRelations[(O,Org)]=Relation[1]

def AddRelations(Orgs1, Orgs2, Relation):
	assert 0
	for O1 in Orgs1:
		for O2 in Orgs2:
				TableRelations[(O1,O2)]=Relation[0]
				TableRelations[(O2,O1)]=Relation[1]
	

def OrganizationsRelatedComplete(Org):	return (OBelowNotCovered[Org],OMaybeCovered[Org],OCovered[Org],OSide[Org],OCovering[Org],OMaybeCovering[Org],OAboveNotCovering[Org])		

def OrganizationsNear(Org,Context=OrganizationsStudied):	return (OMaybeCovered[Org]&Context,OCovered[Org]&Context,OCovering[Org]&Context,OMaybeCovering[Org]&Context)

def OrganizationsRelated(Org):	return (OrganizationsUnder(Org),OSide[Org],OrganizationsAbove(Org))	

