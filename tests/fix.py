from glob import glob
from os import path
from os import rename

old_f = "triclinic"

sys = "tri"

for of in glob(old_f+"/*"):
    fs = of.split(".")
    if "vasp" in of:
        new_f = "symmetry-reduction"
        nf = fs[0].split("/")[1]+"."+fs[1]+"."+sys+"_2_"+fs[2]
    elif "BZM" in of:
        new_f = "BZM"
        nf = fs[0].split("/")[1]+"."+fs[1]+"."+sys+"_"+fs[2]
    else:
        new_f = "symmetry-reduction"
        nf = fs[0].split("/")[1]+"."+fs[1]+"."+sys+"_2_"+fs[2]
    if path.isfile(path.join(new_f, nf)):
        continue
        # print("File exists {}".format(path.join(new_f, nf)))
        # break
    else:
        rename(path.join(of), path.join(new_f, nf))
        
