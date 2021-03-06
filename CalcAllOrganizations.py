# -*- coding: utf-8 -*-
import os
import sys
import gc
import random
import math
import copy


import time

import cPickle
import pickle

import shelve

import itertools


ShowUnion=False#True #False
ShowIntersection=False#True #False
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


def setuparchives(directory):
	global AllStates, AllRndSeeds
	print "I am setting up the ALLSTATES archive"
	#	AllStates=shelve.open("%s%sAllStates.shelve"%(directory,os.sep),writeback=False)     
	#	AllStates=shelve.open("%s%sAllStates1.shelve"%(directory,os.sep))     
	AllStates=shelve.open("%s%sAllStates.shelve"%(directory,os.sep))     
	#if writeback=True the DB gets corrupted when I access all the entries
	AllRndSeeds=shelve.open("%s%sAllRndSeeds.shelve"%(directory,os.sep),writeback=True)
	

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
######## Molecule #########################
###########################################
execfile("%smolecules.py"%chemistry)


###########################################
######## Molecule Library #################
###########################################

def GenerateBasicOrganizations():
	#	 OrgLibrary.check(ReturnAllMolecules())
	OrgLibrary.check(set([])) 


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

def initalisedatabases(DirectoryAC):
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
	reactionDB=MoleculeRactionsDB(DirectoryAC)
	mutationsDB=MoleculeMutationsDB(DirectoryAC)
	LatticeLibrary=LatticeLibraryDB(DirectoryAC)
	GeneratorsLibrary=GeneratorsLibraryDB(DirectoryAC)
	OrganizationName=OrganizationNameDB(DirectoryAC)
	NameOrganization=NameOrganizationDB(DirectoryAC)
	TimeBiggerOrganizationName=TimeBiggerOrganizationNameDB(DirectoryAC)
	IntersectionRelOrg=RelationOrganizationDB("Intersection",DirectoryAC)
	print "I am opening the archives"
	

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

def IsOrganizationSlow(molecules):
	"""
	Returns True if the set of molecules is an organization
	"""
	mol=set(molecules)
	gen=ReactSets(mol,mol)
	if gen==mol:
		return True
	return False

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

def ReactSets(molLeft,molRight):
	"""
	Returns the list of all the molecules generated by a set of molecules, with their relative multiplicity
	"""
	results=set([])
	for m in molLeft:
		for n in molRight:
			res=react(m,n)
			if res!= None:
				results|=set([res])
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
TableIntersection=SymetricTable() #Given two organizations (A,B), stored by name, returns the union of them
TableUnion=SymetricTable() #Given two organizations (A,B), stored by name, returns the union of them. Stores each instance just once unsig the fact that is symmetric

OrganizationsStudied=set([])
#-TableRelations={} #Given two organizations (A,B), returns the relation between them:
##3 =  A greater than B		-3 = B greater than A
##2 = A possibly covers B		-2 = B possibly covers A
##1 = A surely covers B		-1 = B surely covers A
##“I” = Identity (the diagonal)		“<>” = A and B are uncomparable
##			0 = no information


NameSize={}

UnsolvedIntersections=set([])
UnsolvedUnions=set([])

OBelowNotCovered={}
OMaybeCovered={}
OCovered={}
OSide={}
OCovering={}
OMaybeCovering={}
OAboveNotCovering={}

UnionByAUBBuilding=0
UnionTriangulated=0
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


LogOn=False


ToBeSolvedByHandUnionsDict={}
ToBeSolvedByHandIntersectionsDict={}
#Generators={} #the minimum sets that generate something

timeresult=[]
timestart=time.time()
starttime=0.0


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
	

def LogThis(texttolog=None,FileName="Log",directory=".",ForceLog=False):
	if LogOn or ForceLog:
		if texttolog:
			buf=texttolog+"%s"%os.linesep
		else:
			buf=DescribeSituation()+"%s"%os.linesep
		f=open("%s%s%s.txt"%(directory,os.sep,FileName),'a')
		f.write(buf)
		f.close()
	

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

def DescribeSituation():
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

def OrganizationsUnder(Org):				return OCovered[Org]|OMaybeCovered[Org]|OBelowNotCovered[Org]
def OrganizationsAbove(Org):				return OCovering[Org]|OMaybeCovering[Org]|OAboveNotCovering[Org]
def OrganizationsBetween(TopOrg,BottomOrg):	return OrganizationsUnder(TopOrg)&OrganizationsAbove(BottomOrg)
def OrganizationsRelated2(TopOrg,BottomOrg):
	"""Cost:|1|"""
	OrganizationsOverTheTop=OrganizationsAbove(TopOrg)
	OrganizationsUnderTheBottom=OrganizationsUnder(BottomOrg)
	OrganizationsOnSideToBoth=OSide[TopOrg]&OSide[BottomOrg]
	organizationsbetween=OrganizationsBetween(TopOrg,BottomOrg)
	OrganizationsSideTop=OrganizationsAbove(BottomOrg)&OSide[TopOrg] 				#what is on the site the top, but above the bottom
	OrganizationsSideBottom=OSide[BottomOrg]&OrganizationsUnder(TopOrg)
	return (OrganizationsUnderTheBottom,OrganizationsSideBottom,OrganizationsOnSideToBoth,organizationsbetween,OrganizationsSideTop,OrganizationsOverTheTop)
	

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


def GetSetsRelationshipAny(Orgs,R,Context=OrganizationsStudied):
	Result=set([])
	for O in Orgs:
		Result|=R[O]
	return Result&Context

def GetSetsRelationshipAll(Orgs,R,Context=OrganizationsStudied):
	Result=Context
	for O in Orgs:
		Result&=R[O]
	return Result


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


def AddOrganizationOfficially(Org,TopOrg,BottomOrg):
	global OrganizationsStudied
	AddOrganizationRelation(Org,TopOrg,BottomOrg)
	ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections=AddOrganizationUnionIntersection(Org)
	return ToBeSolvedByHandUnions,ToBeSolvedByHandIntersections






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



def SolveByHandUnion(Org):
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


def sortpairs(a,b):
	if a[1]<b[1]:		return -1
	elif a[1]==b[1]:	return 0
	return 1

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


def FindNewOrganizationsAboveByX(BottomOrg,molecules,TopOrg,OrganizationsKnown,SizeSets=2,MaxNumberOrg=None):
	"""Given an organizations BottomOrg, 
	finds all the organizations generated by bottomorg adding one molecule at a time.
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


DownwardMolecules=dict()
UpwardMolecules=dict()
SidewardMolecules=dict()


#from itertools import chain, combinations

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))

def subsets(s):
    return map(set, powerset(s))

	
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





def CheckAllSubsets(Molecules):
	"Checks one by one all the subsets of molecules, to find which one are organisations. But then returns them without storing them"
	organisationsFound=set([])
	subsetsMol=subsets(Molecules)
	for s in subsetsMol:
		print "testing",s
		if IsOrganization(s):
			organisationsFound.add(frozenset(s))
			print "Organisation found", len(organisationsFound), s
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
	
	while OrganisationToExplore in OrganizationsKnown:
		#LastOrgsFound=FindNewOrganizationsAbove(NameBottomOrg,molecules,NameTopOrg)
		LastOrgsFound=FindNewOrganizationsAbove(OrganisationToExplore,molecules,NameTopOrg,OrganizationsKnown,MaxNumberOrg)
		if MaxNumberOrg:
			if len(OrganizationsKnown)>MaxNumberOrg:
				closedatabases()
				sys.exit(0)
		OrganizationsKnown|=LastOrgsFound
		OrganisationToExplore+=1
	#print "MoleculesWithNoUpward=",OrganizationsKnown-set(UpwardMolecules.keys())
	#	for O in UpwardMolecules.keys():	
	#		print "  UpwardMolecules[",O,"]",len(  UpwardMolecules[O]),":",  UpwardMolecules[O]
	
	OrganisationToDeeplyExplore=2
	while OrganisationToDeeplyExplore in OrganizationsKnown:
		#LastOrgsFound=FindNewOrganizationsAbove(NameBottomOrg,molecules,NameTopOrg)
		LastOrgsFound=FindNewOrganizationsAboveByX(OrganisationToDeeplyExplore,molecules,NameTopOrg,OrganizationsKnown,2,MaxNumberOrg)
		if MaxNumberOrg:
			if len(OrganizationsKnown)>MaxNumberOrg:
				closedatabases()
				sys.exit(0)
		OrganizationsKnown|=LastOrgsFound
		OrganisationToDeeplyExplore+=1
	
	
	
	print "OrganizationsKnown=",
	for OrgThis in OrganizationsKnown:
		print "Organisation ",OrgThis,"=",NameOrganization.check(OrgThis)
	
	
	
	#OrgInBetween=FindNewOrganizationsAbove(NameBottomOrg,molecules,NameTopOrg)
	
	
	
	
	#print OrgInBetween






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

