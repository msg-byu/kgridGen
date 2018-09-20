from glob import glob
from os import rename

for fo in glob("unreduced_*"):
    fn = fo.split("unreduced_")[-1]
    rename(fo, fn)
    
