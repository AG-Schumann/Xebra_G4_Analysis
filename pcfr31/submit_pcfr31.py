
#!/usr/bin/python

import sys
import os
from datetime import date

##### INPUT PARAMETER #####
#material_array = []
# Full list of materials: ["LXe", "GXe", "PMT", "Cu", "Ti", "PTFE", "SS", "CryoMat"]
isotope_array = ["Cs137"]
# Full list of Isotopes:  ["U238Pb206", "Co60", "K40", "Cs137", "Ag110m", "Th232Pb208", "U235Pb207", "Ti44Sc44Ca44"]
EVENT_COUNT = 1000000
DATE_STRING = str(date.today())

##### ##### #####

for ISOTOPE_STRING in isotope_array:

    PATH ="/sc/userdata/abismark/output"
    os.makedirs(PATH, exist_ok=True)

    MACRONAME = PATH + "/" + "src_Pointsources_DP" + ".mac"

    f = open(MACRONAME, "w")

    f.write("#VERBOSITY" +'\n' +"/control/verbose 0" +'\n' + "/run/verbose 0" +'\n' +"/event/verbose 0" +'\n' +"/tracking/verbose 0" +'\n' + "/xebra/gun/verbose 0" +'\n' +'\n')
    f.write("#SEED" +'\n' "/run/random/setRandomSeed 0" +'\n' +'\n')
    f.write("# General source settings"  +'\n' +"/xebra/gun/angtype  iso" +'\n' +"/xebra/gun/type Point" +'\n' + "/xebra/gun/center 0. -122. -35.75 mm" +'\n' +"/xebra/gun/particle ion" +'\n' +"/xebra/gun/energy 0 keV" +'\n'+'\n')

    if ISOTOPE_STRING == "Cs137": {f.write("### Cs137" +'\n' +"/xebra/gun/ion 55 137 0 0" +'\n'+'\n')}

    f.write("#ADVANCED RUN OPTIONS" +'\n' + "/Xe/detector/setLXeScintillation false" +'\n' + "/run/writeEmpty false")

    f.close()

    os.system("sbatch -o %s/job_%%j.out simple_job_pcfr31.sh %s %s %i"%(PATH, PATH, ISOTOPE_STRING, EVENT_COUNT))


