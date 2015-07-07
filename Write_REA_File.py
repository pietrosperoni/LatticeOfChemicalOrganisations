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

print"reading lib file: CalcAllOrganizations.py, oops I mean LibraryAC.py"
execfile("LibraryAC.py")


DirectoryAC="ArtificialChemistries%s%s%sExperiments%s%s%s%s"%(os.sep,chemistry,os.sep,os.sep,experiment,os.sep,urllib)
if not os.path.isdir(DirectoryAC):
	os.mkdir(DirectoryAC)

print "data stored in ", DirectoryAC
newallstates=shelve.open("%s%snewallstates.shelve"%(DirectoryAC,os.sep),writeback=False)
print "opening the shelve newallstates"

#setuparchives(directory)
setuparchives(DirectoryAC)
print "setting up the archives"

#DirectoryAC="ArtificialChemistries%s%s%sExperiments%s%s"%(os.sep,chemistry,os.sep,os.sep,experiment)


initalisedatabases(DirectoryAC)
print"initialising the database"



if not os.path.isdir(DirectoryAC):
	os.mkdir(DirectoryAC)
print "data stored in ", DirectoryAC
molecules=ReturnAllMolecules()
InitialiseReaFile(molecules)
print "REA file started"
for m1 in molecules:
	for m2 in molecules:
		m3=react(m1,m2)
print "REA file ended"
closedatabases()
sys.exit(0)