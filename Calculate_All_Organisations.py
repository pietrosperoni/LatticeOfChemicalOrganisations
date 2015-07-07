import time
import os
print"starting"

directory="."
global DirectoryAC
global urllib

execfile("chemistry.py")
print"chemistry=",chemistry
execfile("%s%sArtificialChemistries%s%s%sExperiments%s%s%soptions.py"%(directory,os.sep,os.sep,chemistry,os.sep,os.sep,experiment,os.sep))
DirectoryAC="ArtificialChemistries%s%s%sExperiments%s%s%s%s"%(os.sep,chemistry,os.sep,os.sep,experiment,os.sep,urllib)
if not os.path.isdir(DirectoryAC):
	os.mkdir(DirectoryAC)



execfile("%s%sArtificialChemistries%s%s%smolecules.py"%(directory,os.sep,os.sep,chemistry,os.sep))

print "AC %s loaded"%(urllib)


#DirectoryAC="%s%s%s%s%s"%(directory,os.sep,os.pardir,os.sep,urllib)

#execfile("%s%sArtificialChemistries%s%s%sExperiments%s%s%soptions.py"%(directory,os.sep,os.sep,chemistry,os.sep,os.sep,experiment,os.sep))

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

MapLatticeFromBelow(NameTop, NameBottom)


closedatabases()

sys.exit(0)


