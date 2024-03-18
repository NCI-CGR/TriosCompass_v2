## New improvement for Strelka DNM calls
We had explored to use Strelka to call DNMs, together with DeepVariant and HaplotypeCaller, as demonstrated in the file [Snakefile](./Snakefile).
+ [New changes in the Snakefile](https://github.com/NCI-CGR/TriosCompass_v2/commit/36a6726943af9f0473edfc3150e1121bc2c994ee#diff-47959dfd378b3cd1d39b5515418ee8e4444ab7a6036d5197a5bea82814f928a3):
  + Restrict to “pass” only.
  + Use percentage in filters to address those variants with lower depth.
  + Add new IGV snapshots.
    + New perl script [prepare_yml_from_vcf.pl](./scripts/prepare_yml_from_vcf.pl) to prepare YAML file required by igv_snapshot_maker.
    + Install igv-snapshot-maker to run the updated workflow.
      + pip install -i https://test.pypi.org/simple/ igv-snapshot-maker==1.0.0
+ Updated [summary Excel file](https://github.com/NCI-CGR/TriosCompass_v2/blob/main/cgr_summary.xlsx).
+ New output files under output/call_ism/strelka
  + Users may copy the files (output/call_ism/strelka) to local and start with the Excel files, for example, like output/call_ism/strelka/t0007c1.xlsx
  + The updated DNM calls from Strelka and JIGV should remain in the same location as before under /data/DCEG_Trios/new_cgr_data/TriosCompass_v2.
```bash
### The overall fold struture
tree -L 4  output/call_ism/strelka/t0007c1.xlsx
output/call_ism/
└── strelka
    ├── IGV_Snapshots
    │   ├── t0007c1
    │   │   ├── t0007c1.bat
    │   │   ├── t0007c1_MIE_00001.bat
    │   │   ├── t0007c1_MIE_00001.png
    │   │   ├── t0007c1_MIE_00002.bat
    │   │   ├── t0007c1_MIE_00002.png
    │   │   ├── t0007c1_MIE_00003.bat
    │   │   ├── t0007c1_MIE_00003.png
    │   │   ├── t0007c1_MIE_00004.bat
    │   │   ├── t0007c1_MIE_00004.png
    │   │   ├── t0007c1_MIE_00005.bat
    │   │   ├── t0007c1_MIE_00005.png
    │   │   ├── t0007c1_MIE_00006.bat
    │   │   ├── t0007c1_MIE_00006.png
    │   │   └── t0007c1_ROIs.bat
    │   ├── t0022c1
    │   │   ├── t0022c1.bat
    │   │   ├── t0022c1_MIE_00001.bat
    │   │   ├── t0022c1_MIE_00001.png
    │   │   ├── t0022c1_MIE_00002.bat
    │   │   ├── t0022c1_MIE_00002.png
    │   │   ├── t0022c1_MIE_00003.bat
    │   │   ├── t0022c1_MIE_00003.png
    │   │   ├── t0022c1_MIE_00004.bat
    │   │   ├── t0022c1_MIE_00004.png
    │   │   ├── t0022c1_MIE_00005.bat
    │   │   ├── t0022c1_MIE_00005.png
    │   │   └── t0022c1_ROIs.bat
    │   ├── t0042c1
    │   │   ├── t0042c1.bat
    │   │   ├── t0042c1_MIE_00001.bat
    │   │   ├── t0042c1_MIE_00001.png
    │   │   ├── t0042c1_MIE_00002.bat
    │   │   ├── t0042c1_MIE_00002.png
    │   │   ├── t0042c1_MIE_00003.bat
    │   │   ├── t0042c1_MIE_00003.png
    │   │   ├── t0042c1_MIE_00004.bat
    │   │   ├── t0042c1_MIE_00004.png
    │   │   ├── t0042c1_MIE_00005.bat
    │   │   ├── t0042c1_MIE_00005.png
    │   │   ├── t0042c1_MIE_00006.bat
    │   │   ├── t0042c1_MIE_00006.png
    │   │   ├── t0042c1_MIE_00007.bat
    │   │   ├── t0042c1_MIE_00007.png
    │   │   ├── t0042c1_MIE_00008.bat
    │   │   ├── t0042c1_MIE_00008.png
    │   │   ├── t0042c1_MIE_00009.bat
    │   │   ├── t0042c1_MIE_00009.png
    │   │   ├── t0042c1_MIE_00010.bat
    │   │   ├── t0042c1_MIE_00010.png
    │   │   ├── t0042c1_MIE_00011.bat
    │   │   ├── t0042c1_MIE_00011.png
    │   │   ├── t0042c1_MIE_00012.bat
    │   │   ├── t0042c1_MIE_00012.png
    │   │   └── t0042c1_ROIs.bat
    │   ├── t0058c1
    │   │   ├── t0058c1.bat
...
    ├── t0007c1.xlsx
    ├── t0007c1.yaml
    ├── t0022c1.xlsx
...
    └── t0766c2.yaml

41 directories, 938 files


### Example of the command to transfer data to local (run it at your laptop)
scp -r helix:/data/DCEG_Trios/new_cgr_data/TriosCompass_v2/output/call_ism/strelka ~/Downloads/strelka

### Then open ~/Downloads/strelka/t0007c1.xlsx to review the snapshots
```

---

## Output of TriosCompass with Strelka
+ /data/DCEG_Trios/new_cgr_data/TriosCompass_v2/old_output/old_output

### Summary report
In our latest trios analysis, 107 CGR samples from 40 trios were processed. We developed a Perl script [generate_summary_report.pl](./scripts/generate_summary_report.pl) to generate summary table in Xlxs format.

```bash
### remove the problematic trio t0588c1
ls new_cgr_pedfiles/*.ped  | grep -v t0588c1 > ped.lst

wc -l ped.lst
39 ped.lst

scripts/generate_summary_report.pl cgr_summary.xlsx ped.lst
```

+ [Summary table](./cgr_summary.xlsx)


| FamilyID | ID      | SampleID | FatherSampleID | MotherSampleID | Gender | DNM_Count_DG | DNM_Count_Strelka | JIGV_DG   | JIGV_Strelka |
|----------|---------|----------|----------------|----------------|--------|--------------|-------------------|-----------|--------------|
| t0042    | t0042c1 | SC253873 | SC253877       | SC253872       | M      | 97           | 12                | JIGV HTML | JIGV HTML    |
| t0666    | t0666c1 | SC502245 | SC502256       | SC502234       | F      | 54           | 8                 | JIGV HTML | JIGV HTML    |
| t0705    | t0705c1 | SC742298 | SC742301       | SC742299       | F      | 101          | 8                 | JIGV HTML | JIGV HTML    |
| t0750    | t0750c1 | SC736755 | SC736760       | SC736736       | F      | 93           | 14                | JIGV HTML | JIGV HTML    |
| t0594    | t0594c1 | SC736817 | SC736753       | SC736754       | F      | 110          | 9                 | JIGV HTML | JIGV HTML    |
| t0592    | t0592c1 | SC109501 | SC109499       | SC109498       | F      | 62           | 12                | JIGV HTML | JIGV HTML    |
| t0693    | t0693c2 | SC742220 | SC742221       | SC742222       | F      | 90           | 17                | JIGV HTML | JIGV HTML    |
| t0315    | t0315c2 | SC260715 | SC260727       | SC260729       | M      | 80           | 8                 | JIGV HTML | JIGV HTML    |
| t0575    | t0575c1 | SC109390 | SC109395       | SC109405       | M      | 67           | 12                | JIGV HTML | JIGV HTML    |
| t0766    | t0766c2 | SC736756 | SC736703       | SC730933       | F      | 113          | 11                | JIGV HTML | JIGV HTML    |
| t0271    | t0271c1 | SC109417 | SC109426       | SC109450       | M      | 144          | 11                | JIGV HTML | JIGV HTML    |
| t0600    | t0600c1 | SC109494 | SC109514       | SC109500       | M      | 58           | 6                 | JIGV HTML | JIGV HTML    |
| t0140    | t0140c1 | SC742196 | SC742286       | SC742197       | F      | 88           | 11                | JIGV HTML | JIGV HTML    |
| t0712    | t0712c1 | SC742313 | SC742314       | SC742315       | F      | 114          | 10                | JIGV HTML | JIGV HTML    |
| t0565    | t0565c1 | SC109438 | SC109373       | SC109368       | F      | 68           | 7                 | JIGV HTML | JIGV HTML    |
| t0703    | t0703c1 | SC742295 | SC742296       | SC742297       | F      | 77           | 9                 | JIGV HTML | JIGV HTML    |
| t0483    | t0483c2 | SC742320 | SC742219       | SC742321       | F      | 125          | 8                 | JIGV HTML | JIGV HTML    |
| t0022    | t0022c1 | SC742290 | SC742217       | SC742218       | M      | 64           | 5                 | JIGV HTML | JIGV HTML    |
| t0765    | t0765c1 | SC736700 | SC736738       | SC736739       | M      | 115          | 13                | JIGV HTML | JIGV HTML    |
| t0766    | t0766c1 | SC742179 | SC736703       | SC730933       | M      | 82           | 8                 | JIGV HTML | JIGV HTML    |
| t0739    | t0739c1 | SC736720 | SC736795       | SC736726       | M      | 64           | 11                | JIGV HTML | JIGV HTML    |
| t0575    | t0575c2 | SC109406 | SC109395       | SC109405       | M      | 60           | 11                | JIGV HTML | JIGV HTML    |
| t0679    | t0679c1 | SC742275 | SC742188       | SC742277       | M      | 76           | 9                 | JIGV HTML | JIGV HTML    |
| t0617    | t0617c1 | SC260701 | SC260705       | SC260697       | M      | 131          | 25                | JIGV HTML | JIGV HTML    |
| t0599    | t0599c1 | SC736788 | SC736742       | SC736743       | M      | 106          | 12                | JIGV HTML | JIGV HTML    |
| t0679    | t0679c2 | SC742276 | SC742188       | SC742277       | F      | 105          | 7                 | JIGV HTML | JIGV HTML    |
| t0058    | t0058c2 | SC108472 | SC109353       | SC109336       | M      | 86           | 25                | JIGV HTML | JIGV HTML    |
| t0007    | t0007c1 | SC499427 | SC499423       | SC499428       | M      | 67           | 6                 | JIGV HTML | JIGV HTML    |
| t0058    | t0058c1 | SC109341 | SC109353       | SC109336       | F      | 73           | 16                | JIGV HTML | JIGV HTML    |
| t0565    | t0565c2 | SC109437 | SC109373       | SC109368       | F      | 79           | 8                 | JIGV HTML | JIGV HTML    |
| t0749    | t0749c1 | SC736820 | SC736759       | SC736701       | F      | 90           | 13                | JIGV HTML | JIGV HTML    |
| t0315    | t0315c1 | SC260714 | SC260727       | SC260729       | F      | 74           | 6                 | JIGV HTML | JIGV HTML    |
| t0280    | t0280c1 | SC109512 | SC109516       | SC109507       | M      | 89           | 11                | JIGV HTML | JIGV HTML    |
| t0600    | t0600c2 | SC109495 | SC109514       | SC109500       | F      | 93           | 13                | JIGV HTML | JIGV HTML    |
| t0209    | t0209c1 | SC742216 | SC742271       | SC742272       | F      | 87           | 10                | JIGV HTML | JIGV HTML    |
| t0243    | t0243c1 | SC736787 | SC736702       | SC736772       | F      | 99           | 8                 | JIGV HTML | JIGV HTML    |
| t0707    | t0707c1 | SC742305 | SC742306       | SC742307       | F      | 61           | 8                 | JIGV HTML | JIGV HTML    |
| t0760    | t0760c1 | SC736729 | SC736824       | SC736803       | F      | 106          | 9                 | JIGV HTML | JIGV HTML    |
| t0623    | t0623c1 | SC253881 | SC253893       | SC253894       | F      | 52           | 6                 | JIGV HTML | JIGV HTML    |


:bookmark: In the Excel file, the last two columns of "JIGV HTML" are built in with HTTP links to the local HTML pages under old_output/call_JIGV.  Therefore, the HTTP links will not work properly unless the report file is located in the right position. 

```bash
cgr_summary.xlsx
old_output/call_JIGV/
├── D_and_G_t0007c1.JIGV.html
├── D_and_G_t0022c1.JIGV.html
...
├── D_and_G_t0766c1.JIGV.html
├── D_and_G_t0766c2.JIGV.html
├── strelka_t0007c1.JIGV.html
├── strelka_t0022c1.JIGV.html
...
├── strelka_t0766c1.JIGV.html
└── strelka_t0766c2.JIGV.html
```

---
