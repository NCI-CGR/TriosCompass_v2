# Changelog

### [1.0.1](https://www.github.com/NCI-CGR/TriosCompass_v2/compare/v1.0.0...v1.0.1) (2026-05-18)


### Bug Fixes

* add memory requirements to GATK CombineGVCFs step ([c24b2fe](https://www.github.com/NCI-CGR/TriosCompass_v2/commit/c24b2fea46ee670a9952f79adff2e3c1f22d7360))
* **workflow/profiles/slurm/config.yaml:** change memory for glnexus_dv and aggregate_phase ([55605b1](https://www.github.com/NCI-CGR/TriosCompass_v2/commit/55605b1325f6716618e359d23df4dbba46068c50))
* **workflow/profiles/slurm/config.yaml:** for the slurm setting used in CCAD (some of them might not be compatible at Biowulf) ([b6733d9](https://www.github.com/NCI-CGR/TriosCompass_v2/commit/b6733d99f86aa2ec938ffc542bda033950472322))
* **workflow/report/dnm.rst:** fix the typo in rst ([957f593](https://www.github.com/NCI-CGR/TriosCompass_v2/commit/957f593de2078334a088d51314b946f489db4956))
* **workflow/rules/dnSTR.smk:** dnSTR summary for each trio (not family) ([8bc1426](https://www.github.com/NCI-CGR/TriosCompass_v2/commit/8bc1426e3d43bb12172aec7f767291296c39ff45))
* **workflow/rules/dnSV_manta.smk:** replace manta docker image to fix: taskExitCode -11 ([0513d22](https://www.github.com/NCI-CGR/TriosCompass_v2/commit/0513d22322b1ec93d46708dbb6a326c42c330bd9))
* **workflow/rules/hipstr.smk:** exit due to zgrep and also fix one typo ([2d1746e](https://www.github.com/NCI-CGR/TriosCompass_v2/commit/2d1746e427acf6de32472e754b65293358d3e5e8))
* **workflow/rules/multiqc.smk:** use docker://multiqc/multiqc:v1.33 to fix issue in the use of multiqc conda or wrapper ([6b5659d](https://www.github.com/NCI-CGR/TriosCompass_v2/commit/6b5659d6c3bba521afa5c45e323934bdce4706da))
* **workflow/rules/pedigree.smk:** sorted pedigree file in the order: father,mother,child ([b92463c](https://www.github.com/NCI-CGR/TriosCompass_v2/commit/b92463c09dd9f9007ac520e0db9e3dec6ae735af))
