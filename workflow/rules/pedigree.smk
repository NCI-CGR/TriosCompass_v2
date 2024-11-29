import peds

### Define trios
ped_dir = config["ped_dir"]
ped_files = glob.glob(config["ped_dir"] + "/*.ped")

families = {}
for fn in ped_files:
    f=peds.open_ped(fn)[0]
    families[f.id]=f

fam_ids = list(families.keys()) 

child_ids = [[person.id for person in families[fid] if families[fid].get_father(person) ][
0] for fid in fam_ids]

final_subjs = list(set([p.id for f in fam_ids for p  in families[f] ]))

CHILD_DICT=dict(zip(fam_ids,child_ids))