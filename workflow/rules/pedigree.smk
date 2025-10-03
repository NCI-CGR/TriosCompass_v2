import peds

### Define trios
ped_dir = config["ped_dir"]
ped_files = glob.glob(config["ped_dir"] + "/*.ped")

# unsorted families
_families = {}
for fn in ped_files:
    f=peds.open_ped(fn)[0]
    _families[f.id]=f

### sorted families
# https://github.com/jeremymcrae/peds/blob/master/tests/test_link_members.py
# assume ped is defined for trios only and only children have both parents defined
# sort families to father, mother, child and for convenient dnSV calling
families={}
for family_id, family_obj in _families.items():
    child = [person for person in family_obj if family_obj.get_father(person)][0]

    
    new_fam = peds.Family(family_id)
    for p in [family_obj.get_father(child), family_obj.get_mother(child), child]:
        new_fam.add_person(p)
        
    families[family_id] = peds.ped.link_members(new_fam)

fam_ids = list(families.keys()) 

child_ids = [[person.id for person in families[fid] if families[fid].get_father(person) ][0] for fid in fam_ids]

final_subjs = list(set([p.id for f in fam_ids for p  in families[f] ]))

CHILD_DICT=dict(zip(fam_ids,child_ids))