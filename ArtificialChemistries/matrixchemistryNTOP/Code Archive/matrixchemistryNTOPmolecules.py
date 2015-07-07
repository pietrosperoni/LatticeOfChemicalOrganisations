# -*- coding: utf-8 -*-
#import numpy
##import hashlib; import base64
import binascii
global urllib

global DirectoryAC


def writeoptions(directory): #unused
	"""a function that writes a file options.py with all the value of the options.
	The options are stored as python variables, and it is enough to execute the options.py program to get those options"""
	options=""
	print "writing options"
	global urllib
	global DirectoryAC
	options+="AC=\"matrixchemistry\"%s"%os.linesep
	options+="MainToPrimary=%s%s"%(MainToPrimary,os.linesep)
	options+="MainToSecondary=%s%s"%(MainToSecondary,os.linesep)
	options+="Result2Main=%s%s"%(Result2Main,os.linesep)
	options+="IllegalMolecules=%s%s"%(IllegalMolecules,os.linesep)
	#	m = hashlib.md5().update(options)
	#	urllib=base64.urlsafe_b64encode(m.digest())
	urllib="matrixchemistry_%s"%(binascii.crc32(options))
	libdir="%s%s%s%s%s"%(directory,os.sep,os.pardir,os.sep,urllib)
	#directoryAC=libdir

	if not os.path.isdir(libdir):
			os.mkdir(libdir)
	open("%s%sACType.py"%(libdir,os.sep),'w').write(options)

	options+="urllib='%s'%s"%(urllib,os.linesep)
	options+="mutation_rate_bit_flip=%s%s"%(mutation_rate_bit_flip,os.linesep)
	options+="starting_length=%s%s"%(starting_length,os.linesep)
	options+="starting_cell_size=%s%s"%(starting_cell_size,os.linesep)
	options+="LengthofExperiments=%s%s"%(LengthofExperiments,os.linesep)
	options+="TimeBetweenSavings=%s%s"%(TimeBetweenSavings,os.linesep)
	open("%s%smatrixchemistryoptions.py"%(directory,os.sep),'w').write(options)









#mutationsDB=MoleculeMutationsDB()

def testrightidentity(molecule,context):
	for m in context:
		if reactionDB.check(m,molecule)!=m:
			return False
	return True

def testleftidentity(molecule,context):
	for m in context:
		if reactionDB.check(molecule,m)!=m:
			return False
	return True

def testrightzero(molecule,context):
	for m in context:
		if reactionDB.check(m,molecule)!=molecule:
			return False
	return True

def testleftzero(molecule,context):
	for m in context:
		if reactionDB.check(molecule,m)!=molecule:
			return False
	return True


	

	
###########################################
######## Molecule #########################
###########################################

def name2valuesize(name):
	"calculates and returns the value and the size from the name"
	size=int(math.floor(math.log(name,2)))
	return int(name-2**size),int(size)



def colormolecule1(molname):
	saturation=0.5
	if reactionDB.check(molname,molname)==molname:
		hue=0.0
	else:
		hue=0.66666
	lightness=0.9
	col=colorsys.hsv_to_rgb(hue, saturation, lightness)
	return col

def colormolecule2(molname):
	mol=moleculelibrary.check(molname)
	red=mol.mainlist[0]*128+mol.mainlist[1]*64+mol.mainlist[1]*32
	green=mol.mainlist[3]*128+mol.mainlist[4]*64+mol.mainlist[5]*32
	blue=mol.mainlist[6]*128+mol.mainlist[7]*64+mol.mainlist[8]*32
	
	return (red/256.0,green/256.0,blue/256.0)

def ColorLayout():
	while 1:
		yield colormolecule1
		yield colormolecule2


def writereact(t,h,option,molecules,numbermolecules):
	op=option%6
	if not op:   return ""
	if op==1 :
		res=reactionDB.check(t,h)
		if res:  return "%i\t"%res
		else  :  return "\t"
	if op==2 :
		res=reactionDB.check(t,h)
		if t==h:
			if   res==t: return "*\t"
			elif res   : return "%i\t"%res
			else	   : return "\t"
		else:
			if   res==t: return ">\t"
			elif res==h: return "^\t"
			elif res   : return "%i\t"%res
			else	   : return "\t"
	if op==3 :
		res=reactionDB.check(t,h)
		if t==h:
			if   res==t:		return "*\t"
			elif res   :		return "%i\t"%res
			else	   :		return "\t"
		else:
			ser=reactionDB.check(h,t)
			if   res==t:
				if   ser==t:	return "*->\t"
				elif ser==h:	return "==>\t"
				elif ser:	   return "->?\t"
				else :		  return " ?\t"
			elif   res==h:
				if   ser==t:	return "==^\t"
				elif ser==h:	return "*^\t"
				elif ser:	   return "^?\t"
				else :		  return "^\t"
			elif res:
				if   ser==t:	return "%i^\t"%res
				elif ser==h:	return "%i>\t"%res
				elif ser:
					if ser==res:return "*%i\t"%res
					else:	   return "%i?\t"%res
				else :		  return "%i\t"%res
			else:
				if ser==t:	  return " ^\t"%res
				elif ser==h:	return " >\t"%res
				elif ser:	   return " ?\t"%res
				else:		   return "\t"
	if op==4 :
		res=reactionDB.check(t,h)
		try:						ind=molecules.index(res)+1
		except ValueError:		  ind=0
		if t==h:
			if   res==t:					return "*\t"
			elif res:
				if ind==0:				  return "NEW\t"
				elif ind<=numbermolecules:   return "# %i\t"%ind
				else:					   return "%i\t"%res
			else:						   return "\t"
		else:
			if   res==t:					return ">\t"
			elif res==h:					return "^\t"
			elif res:
				if ind==0:				  return "NEW\t"
				elif ind<=numbermolecules:   return "# %i\t"%ind
				else:					   return "%i\t"%res
			else:						   return "\t"
	if op==5 :
		res=reactionDB.check(t,h)
		try:						ind=molecules.index(res)+1
		except ValueError:		  ind=0
		if t==h:
			if res:
				if ind==0:				  return "NEW\t"
				elif ind<=numbermolecules:   return "# %i\t"%ind
				else:					   return "%i\t"%res
			else:						   return "\t"
		else:
			if res:
				if ind==0:				  return "NEW\t"
				elif ind<=numbermolecules:   return "# %i\t"%ind
				else:					   return "%i\t"%res
			else:						   return "\t"



	   
## moved to aclib making it general
##def GeneratedMolecules(molecules):
##	"""
##	Returns the list of all the molecules generated by a set of molecules, with their relative multiplicity
##	"""
##	results=mydict()
##	for m in molecules:
##		for n in molecules:
##			res=reactionDB.check(m,n)
##			if res!= None:
##				try:
##					results[res]+=1
##				except KeyError:
##					results[res]=1
##	return results
##
##


