# Changelog

## [0.2.0](https://github.com/RIVM-bioinformatics/juno-mapping/compare/v0.1.0...v0.2.0) (2023-12-18)


### Features

* test data and unit tests ([234650f](https://github.com/RIVM-bioinformatics/juno-mapping/commit/234650f70ec75f8745ffae59ce5db0b64509f043))


### Bug Fixes

* correct syntax when using --disable-mask ([9e1834f](https://github.com/RIVM-bioinformatics/juno-mapping/commit/9e1834ffa768916d6f0973d9cdff3459ce4731ee))

## 0.1.0 (2023-07-21)


### Features

* add filter status of variants to mqc ([2140f48](https://github.com/RIVM-bioinformatics/juno-mapping/commit/2140f48122592b8ee8b52c1903976c1b203c4a80))
* add mapping and variant calling ([2f46ea8](https://github.com/RIVM-bioinformatics/juno-mapping/commit/2f46ea8606305a44729458159255574f347cc344))
* add multiqc ([6e42991](https://github.com/RIVM-bioinformatics/juno-mapping/commit/6e4299175114d082ae0c67eb61109264cb5004d6))
* add multiqc report ([402b7f0](https://github.com/RIVM-bioinformatics/juno-mapping/commit/402b7f07f2f55b495011cb82b8128c000ad53d13))
* add numeric argument validation ([efdb13b](https://github.com/RIVM-bioinformatics/juno-mapping/commit/efdb13bd5f5082147826dc1cc85581866cfbf660))
* add options to wrapper ([14ae8a5](https://github.com/RIVM-bioinformatics/juno-mapping/commit/14ae8a561060ad0fd8fa1ab6f8bbc4ca7978fb7c))
* add qc measures ([90c4248](https://github.com/RIVM-bioinformatics/juno-mapping/commit/90c42482ef57b35c5b761d818ce9baad85f37a1b))
* add resource management for hpc ([d498478](https://github.com/RIVM-bioinformatics/juno-mapping/commit/d49847845c35425bcd7ff292c01e4c0010b252c9))
* add workflow shared with juno-assembly ([8a3246e](https://github.com/RIVM-bioinformatics/juno-mapping/commit/8a3246e34ad8bc21a014594ac204753cefb304af))
* allow custom reference ([fa0d89d](https://github.com/RIVM-bioinformatics/juno-mapping/commit/fa0d89d0ec26f8ed1f4df94247a4fc3b067cd5a8))
* wrapper script ([c75d4fb](https://github.com/RIVM-bioinformatics/juno-mapping/commit/c75d4fb25f986fb976e30b2dfc9653a48fbb72ab))


### Bug Fixes

* add another auxiliary script ([c79a68a](https://github.com/RIVM-bioinformatics/juno-mapping/commit/c79a68aeda269a9e882df6928501c53a70f4c0a0))
* add auxiliary scripts ([0db0a10](https://github.com/RIVM-bioinformatics/juno-mapping/commit/0db0a1031bbe6d512d5dde7e6b45f2c2f2ed3d0b))
* allow disabling of masking without errors ([6dc7576](https://github.com/RIVM-bioinformatics/juno-mapping/commit/6dc75762f5bbf0939f1ea45e301577c530362d43))
* class name in main function ([5c5a402](https://github.com/RIVM-bioinformatics/juno-mapping/commit/5c5a4026a4361eceea25e265102cd46b820c3f5b))
* conda env path kraken2 ([874109a](https://github.com/RIVM-bioinformatics/juno-mapping/commit/874109a703c8cc6eae58595e1a65f74eebb24ecb))
* conda envs for all rules ([77cb6e4](https://github.com/RIVM-bioinformatics/juno-mapping/commit/77cb6e427c882eb7ab3f3a3297cbf8a9b4ec8b75))
* copy mask file to output ([c9a5ad6](https://github.com/RIVM-bioinformatics/juno-mapping/commit/c9a5ad646a10e17cc27a7e8b351fbd49135c9a8d))
* default variant depth ([51a4db2](https://github.com/RIVM-bioinformatics/juno-mapping/commit/51a4db28fc2780806fa43c530d879c8d86a03e8f))
* gitignore not ignoring config files ([748f19b](https://github.com/RIVM-bioinformatics/juno-mapping/commit/748f19bc3e07067d2b0362821bf8999b675ed5dc))
* import argparse earlier because of typing ([f33fd7c](https://github.com/RIVM-bioinformatics/juno-mapping/commit/f33fd7ce97a72ca81db7ecb0711132547039fd72))
* keep k2report ([8954c3e](https://github.com/RIVM-bioinformatics/juno-mapping/commit/8954c3e2d18688c2f914efc4bc0b936b42b5ce4d))
* log paths ([7dd870f](https://github.com/RIVM-bioinformatics/juno-mapping/commit/7dd870ff0f9deeee457c7a0a46385068aff5e74a))
* min_af parsing ([f909e59](https://github.com/RIVM-bioinformatics/juno-mapping/commit/f909e59075419de59a2375773a016fc9550d5fa9))
* parsing of disable-mask config value ([fbda34b](https://github.com/RIVM-bioinformatics/juno-mapping/commit/fbda34bc2f76a0bcbf33e9e3bbec96eb82b07450))
* reorder filtering of variants ([645fd6e](https://github.com/RIVM-bioinformatics/juno-mapping/commit/645fd6e524582b5166f14425a087dd62cf38b247))
* resources for species summary rule ([152d57e](https://github.com/RIVM-bioinformatics/juno-mapping/commit/152d57e5fcdea46f68f72a095686c65d2fb75c7c))
* rework mamba env activation ([b368a64](https://github.com/RIVM-bioinformatics/juno-mapping/commit/b368a64957d2f9a090fd1819b759f77278d6353d))
* rule all file suffix ([38a2fa3](https://github.com/RIVM-bioinformatics/juno-mapping/commit/38a2fa3ca93beff0d08b28e5fae4b98cd397b36e))
* samtools stats rule paths ([0773dd0](https://github.com/RIVM-bioinformatics/juno-mapping/commit/0773dd0354379fa276154b0c91c58e0b3f100bbf))
* test gitignore ([5076a04](https://github.com/RIVM-bioinformatics/juno-mapping/commit/5076a04266cc71ef29d89f3e89c175a31cbd1eb9))
* test gitignore again ([4df2593](https://github.com/RIVM-bioinformatics/juno-mapping/commit/4df2593f970728e3c40edb914c3c86fcd62a4250))
* typo ([ea009ab](https://github.com/RIVM-bioinformatics/juno-mapping/commit/ea009ab8627c6735a12c5ebc9456883686abd294))
* typo in masking ([e9fbaf3](https://github.com/RIVM-bioinformatics/juno-mapping/commit/e9fbaf35689c1f1f5e8b3830b31576e441ce7608))
* typo in select_snps rule ([4045dc4](https://github.com/RIVM-bioinformatics/juno-mapping/commit/4045dc4274153578e75c45969d2b80efe926fa5d))
* working ignore ([95bc2a3](https://github.com/RIVM-bioinformatics/juno-mapping/commit/95bc2a3af4f331fe900ff762be8e53cada1d26c4))


### Dependencies

* add conda defaults ([60c53b5](https://github.com/RIVM-bioinformatics/juno-mapping/commit/60c53b5e5acb45e8ad44a2cea5a4ac2e441ae46d))
* add containers for all rules ([cb028d1](https://github.com/RIVM-bioinformatics/juno-mapping/commit/cb028d1ec10308c7e25c21a3a5585c2c573fe482))
* remove anaconda and defaults channel and add no defaults channel ([f6d2f05](https://github.com/RIVM-bioinformatics/juno-mapping/commit/f6d2f05bf94e853fca975d178347783ea9948630))


### Documentation

* update README ([51aa80c](https://github.com/RIVM-bioinformatics/juno-mapping/commit/51aa80c52b18e9e07a87369be6680a07c7553364))
* update version.py ([14e3edf](https://github.com/RIVM-bioinformatics/juno-mapping/commit/14e3edf334749377a7906995d5bfa0688ed9baad))
