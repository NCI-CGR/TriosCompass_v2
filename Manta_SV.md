# Call dnSV with *manta*

## Introduction
In this revision, dnSV callign by *manta* is incorporated to TriosCompass.  The original dnSV calling approach (i.e., GraphTyper2+smoove) is still retained for the time being.  Users may select to enable/disable dnSV calling using [the configure yaml file](./config/GIAB_40X_dnSV.yaml):

```yaml
dnSV:
  enable: True
  GraphTyper2+smoove:
    enable: True
    exclude_bed: "ref/exclude.cnvnator_100bp.GRCh38.20170403.bed"
  Manta:
    enable: True

```

The output of manta dnSV will be available at:
+ {output_dir}/
  + dnSV_manta/ 
  + dnSV_manta_summary/

The launch script is similar as before except a new command-line argument "--latency-wait 600", which is to warrant sufficient latency time for *manta* to generate the targeted output.