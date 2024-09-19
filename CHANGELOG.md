# Changelog

## [0.4.5](https://github.com/RIVM-bioinformatics/juno-mapping/compare/v0.4.4...v0.4.5) (2024-09-19)


### Bug Fixes

* correct metadata var ([0cc6f3d](https://github.com/RIVM-bioinformatics/juno-mapping/commit/0cc6f3d11bcbc9692df4e51dac26c41bd4e439c6))
* propagate user uploaded metadata ([0fe1218](https://github.com/RIVM-bioinformatics/juno-mapping/commit/0fe12185606e3b81ee1a6da7601b4d4c1fc945c8))
* quote pattern ([876945f](https://github.com/RIVM-bioinformatics/juno-mapping/commit/876945fbffd55b32d847f4e7db6ab3377ad46dc0))

## [0.4.4](https://github.com/RIVM-bioinformatics/juno-mapping/compare/v0.4.3...v0.4.4) (2024-09-12)


### Bug Fixes

* prevent unbound var error on irods acc ([9c89b0f](https://github.com/RIVM-bioinformatics/juno-mapping/commit/9c89b0fbf9e6dd3a6c5174f511ace6413bd5c584))

## [0.4.3](https://github.com/RIVM-bioinformatics/juno-mapping/compare/v0.4.2...v0.4.3) (2024-07-25)


### Bug Fixes

* update cli ([5304d85](https://github.com/RIVM-bioinformatics/juno-mapping/commit/5304d85313523a883966e93344854a2ee45398f2))


### Dependencies

* update juno library ([710bf39](https://github.com/RIVM-bioinformatics/juno-mapping/commit/710bf39d9817157cafddba5a5e32ce83cc294335))

## [0.4.2](https://github.com/RIVM-bioinformatics/juno-mapping/compare/v0.4.1...v0.4.2) (2024-07-01)


### Bug Fixes

* increase mqc mem ([6a0ce31](https://github.com/RIVM-bioinformatics/juno-mapping/commit/6a0ce31bf095daec364ec9de2a5171cf8ecbb569))
* split variants into biallelics ([be3c319](https://github.com/RIVM-bioinformatics/juno-mapping/commit/be3c3191b8bc1940cdc68b14afcaf1d9c40ce812))

## [0.4.1](https://github.com/RIVM-bioinformatics/juno-mapping/compare/v0.4.0...v0.4.1) (2024-05-14)


### Bug Fixes

* improve repeatability ([0f9e3eb](https://github.com/RIVM-bioinformatics/juno-mapping/commit/0f9e3eb02c2232a4da1995dddf12f9c0b3a7d724))

## [0.4.0](https://github.com/RIVM-bioinformatics/juno-mapping/compare/v0.3.2...v0.4.0) (2024-04-16)


### Features

* call large deletions using delly ([bd06022](https://github.com/RIVM-bioinformatics/juno-mapping/commit/bd060229e5b13d05a3a11556d048315449dccae1))


### Bug Fixes

* rm unneeded check ([6c30463](https://github.com/RIVM-bioinformatics/juno-mapping/commit/6c30463ca9607f804f1cde9a087f00719ea8072b))

## [0.3.2](https://github.com/RIVM-bioinformatics/juno-mapping/compare/v0.3.1...v0.3.2) (2024-03-27)


### Bug Fixes

* default mtb ref path ([a899e93](https://github.com/RIVM-bioinformatics/juno-mapping/commit/a899e93d5fcc7fdc61f891b6c5595d9872d866bc))
* set hpc resource config ([a2e64a3](https://github.com/RIVM-bioinformatics/juno-mapping/commit/a2e64a341dc3de220ee884e75e494d5b9a699678))

## [0.3.1](https://github.com/RIVM-bioinformatics/juno-mapping/compare/v0.3.0...v0.3.1) (2024-02-07)


### Bug Fixes

* distinguish soft and hard AF filter ([ee23734](https://github.com/RIVM-bioinformatics/juno-mapping/commit/ee2373473fec6a4d510d1fa323aec21a560e3acd))
* fix paths when not masking regions ([8e8a850](https://github.com/RIVM-bioinformatics/juno-mapping/commit/8e8a8500635db7608a6f9d585cb3ba3d13db3e99))

## [0.3.0](https://github.com/RIVM-bioinformatics/juno-mapping/compare/v0.2.1...v0.3.0) (2024-01-25)


### Features

* add filtering of variants based on strand bias ([3be232a](https://github.com/RIVM-bioinformatics/juno-mapping/commit/3be232abd10f27d3052a81c20cf89f9de8188564))


### Dependencies

* smk to 7.32.0 ([4a59f6a](https://github.com/RIVM-bioinformatics/juno-mapping/commit/4a59f6a6f276b619d16415d0ea2463cce4a61e35))
* update master conda env ([7ebcb10](https://github.com/RIVM-bioinformatics/juno-mapping/commit/7ebcb107dd497a1a1b0d1be5aa79411f7a791600))

## [0.2.1](https://github.com/RIVM-bioinformatics/juno-mapping/compare/v0.2.0...v0.2.1) (2023-12-19)


### Bug Fixes

* conda compatibility of picard ([6803b0c](https://github.com/RIVM-bioinformatics/juno-mapping/commit/6803b0c943320936999ad5fd895284d680c7f4d1))


### Dependencies

* add conda envs ([bb1134a](https://github.com/RIVM-bioinformatics/juno-mapping/commit/bb1134abbdbb527bc606833ab7dbbd989a7bf51a))
* bump multiqc ([94eff56](https://github.com/RIVM-bioinformatics/juno-mapping/commit/94eff564a0e93bfa9b7b78f71827e8dfc157e21a))
* channel reorder ([7564701](https://github.com/RIVM-bioinformatics/juno-mapping/commit/75647010092f41abefd608d759bdc1a530b8c386))
* reorder channels ([293ef1b](https://github.com/RIVM-bioinformatics/juno-mapping/commit/293ef1b155a8df6afe8d71b148a245bbf6b79368))
* reorder channels for openssl bug ([fb8a21f](https://github.com/RIVM-bioinformatics/juno-mapping/commit/fb8a21f941925b1695832ec747a0374186dd54d8))
* set nodefaults for conda envs ([f662649](https://github.com/RIVM-bioinformatics/juno-mapping/commit/f66264991de330f3c702677d3f0cd26d6f5eceab))

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
