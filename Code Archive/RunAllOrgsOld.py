import time
import os
print"starting"

directory="."
global DirectoryAC
global urllib

execfile("chemistry.py")
print"chemistry=",chemistry
execfile("%s%sArtificialChemistries%s%s%sExperiments%s%s%soptions.py"%(directory,os.sep,os.sep,chemistry,os.sep,os.sep,experiment,os.sep))
execfile("%s%sArtificialChemistries%s%s%smolecules.py"%(directory,os.sep,os.sep,chemistry,os.sep))

print "AC %s loaded"%(urllib)

#DirectoryAC="%s%s%s%s%s"%(directory,os.sep,os.pardir,os.sep,urllib)
DirectoryAC="ArtificialChemistries%s%s%sExperiments%s%s%s%s"%(os.sep,chemistry,os.sep,os.sep,experiment,os.sep,urllib)

#execfile("%s%sArtificialChemistries%s%s%sExperiments%s%s%soptions.py"%(directory,os.sep,os.sep,chemistry,os.sep,os.sep,experiment,os.sep))
if not os.path.isdir(DirectoryAC):
	os.mkdir(DirectoryAC)

print "data stored in ", DirectoryAC
#execfile("AClib.py")

print"reading lib file: CalcAllOrganizations.py, oops I mean LibraryAC.py"
#execfile("CalcAllOrganizations.py")
execfile("LibraryAC.py")


if len(sys.argv)>1:
	pass
#	 =int(sys.argv[1]) 
#	 if len(sys.argv)>2:
#		 =int(sys.argv[2]) 

################################################################################
actual=set()
molecules={}
print"initialising the molecule variable"

initalisedatabases(DirectoryAC)
print"initialising the database"

#newallstates=shelve.open("%s%snewallstates.shelve"%(directory,os.sep),writeback=False)
newallstates=shelve.open("%s%snewallstates.shelve"%(DirectoryAC,os.sep),writeback=False)
print "opening the shelve newallstates"

#setuparchives(directory)
setuparchives(DirectoryAC)
print "setting up the archives"

molecules=ReturnAllMolecules()

print "those are the molecules"
print molecules


InitialiseReaFile(molecules)
print "REA file initialised"

TopSet=set(molecules)
TopOrganization=Organization(TopSet)

BottomSet=set([])
BottomOrganization=Organization(BottomSet)

NameTop=OrganizationName.check(TopOrganization)           #OrganizationName is the OrganizationNameDB opened  in initalisedatabases(...)
NameBottom=OrganizationName.check(BottomOrganization)

print "NameTop=",NameTop
print "NameBottom=",NameBottom

print "going to study the Lattice"
sys.stdout.flush()


#brutforce:
#OrganizationsFoundStudyingAll=CheckAllSubsetsBySize(set(molecules))

MapLatticeFromBelow(NameTop, NameBottom)
#MapLattice(NameTop, NameBottom)
#MapLattice(NameTop, NameBottom,5000)

#OrganizationsFoundStudyingAll=CheckAllSubsets(set(molecules))





closedatabases()

sys.exit(0)





#for s in AllStates.keys():
#	newallstates[s]=AllStates[s]

for m in mol.keys():
	temp=numpy.mat(mol[m])
	molecules[m]=temp.T

states=getstates()
#print states

#newallstates={}
#newstepstext=newallstates.keys()
#print "newstepstext" ,newstepstext
#newstepstext.sort()
#newstepsint=[int(d) for d in newstepstext]
#print newstepsint

#for s in newallstates.keys():
#	pass
	#print newallstates[s]

print "number states recorded in new=",len(newallstates)


#for whichstep in newstepstext:
#	print "whichstep",whichstep
#	print newallstates[whichstep]
#	newallstates.sync()

#print "number states recorded in new=",len(newallstates)


#for whichstep in newstepsint:
#	print "whichstep",whichstep
#	print newallstates["%10s"%whichstep]

	
#print newallstates

def readstate (whichstep):
	moleculestemp=AllStates["%10s"%whichstep]
#	b=Cell(AllStates["%10s"%whichstep])
	print "moleculestemp=", moleculestemp
	b=Cell(moleculestemp)
	return b


print "number states recorded=",len(states)
#print "number states recorded in new=",len(newallstates)


#for s in states:
#	b=readstate(s)
#	print "molecules=",b.moleculesdict.keys()
#	if(not IsOrganization(b.moleculesdict.keys())):
#		print "is not an organization, the organization is:"
#		orgpresent=OrgLibrary.check(b.moleculesdict.keys())
#		print orgpresent
#		print "organization known:", GetNKnownOrg()

#AllOrg=set(NameOrganization.DB.keys())
AllOrg=NameOrganization.DB.keys()
NOrg=len(AllOrg)
print NOrg
Buf=""
#for Pointer in range(NOrg):
for Pointer in range(600):
#	print "hi"
	print Pointer,
#	print "hi2"
	O1=AllOrg[Pointer]
#	print "hi3"
#	print AllOrg[Pointer]
#	print "hi4"	
	print len(NameOrganization.DB.keys())	
#	print "hi5"	
	for SecondPointer in range(Pointer+1,NOrg):#		RelOrg1Org2(A,B)
#		Buf+="%s "%IsOrg1inOrg2(O1,AllOrg[SecondPointer])
#		Buf+="%s "%RelOrg1Org2(O1,AllOrg[SecondPointer])
		#Buf+="%s "%Org1intersectionOrg2(O1,AllOrg[SecondPointer])	
		Org1intersectionOrg2(O1,AllOrg[SecondPointer])	
	if Buf:
		print Buf
		Buf=""

print "hi6"	
	
print "NOrg=",len(NameOrganization.DB.keys())

print "hi7"	

#for O1O2 in findsubsets(AllOrg,2):
#	print O1O2,
#for Organiz1 in AllOrg:
#	for Organiz2 in AllOrg:
#		print IsOrg1inOrg2(Organiz1,Organiz2),



#WriteLatticesKnownInGeneral()



#indexstep=0
#tinymove=0
#print "pippo"

#print "all states=",states

#fromtimestep=states[0]

#totimestep=states[-1]

#stepstudied=0
#diversity=[getdiversity(s) for s in states]
print "hi8"	

diversity=getdiversitylist()
print "hi9"	

#print AllStates
closedatabases()
#print len(AllStates)
print 'hi10 \a'	
