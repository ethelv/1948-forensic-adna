#!/usr/bin/env python
"""
Select the reference individuals kept for ADMIXTURE.

The full AADR reference panel is dominated by European populations, which would
otherwise absorb most ADMIXTURE components. We drop the European / Caucasus
populations and keep the Middle-Eastern, North-African and Jewish groups, then
write the individual list used to subset the merged EIGENSTRAT/PLINK dataset.

Input (not shipped): merged.ind (target + AADR reference, EIGENSTRAT).
Output: iids_to_keep (one individual ID per line) -> feed to `plink --keep`.
"""
import pandas as pd

ind = pd.read_csv("merged.ind", delim_whitespace=True, names=["name", "sex", "pop"])

# European / Caucasus populations removed for the ADMIXTURE analysis
to_remove = {
    'Basque.HO','Scottish.HO','Sicilian.HO','Spanish.HO','Spanish_North.HO','Icelandic.HO',
    'Ukrainian.HO','Romanian.HO','Russian.HO','English.HO','Estonian.HO','Finnish.HO','French.HO',
    'Belarusian.HO','Bulgarian.HO','Mordovian.HO','Sardinian.HO','Adygei.HO','Czech.HO','Croatian.HO',
    'Lithuanian.HO','Italian_North.HO','Hungarian.HO','Norwegian.HO','Balkar.HO','Georgian.HO',
    'Orcadian.HO','Chechen.HO','Ossetian.HO','Cypriot.HO','Lezgin.HO','Albanian.HO','Abkhasian.HO',
    'Armenian.HO','Assyrian.HO','Kumyk.HO',
}
keep = ind[~ind["pop"].isin(to_remove)]
keep[["name"]].to_csv("iids_to_keep", index=False, header=False, sep="\t")
print(f"kept {len(keep)} individuals in {keep['pop'].nunique()} populations -> iids_to_keep")
