# Workflow testing

This folder contains scripts and expected results for testing the four potential workflow combinations. These are:

* Illumina + Viridian workflow
* Illumina + Artic workflow
* ONT + Viridian workflow
* ONT + Artic workflow

## Running the tests

Run the `GPAS_tests.bash` script from the headnode after SSH login.

```bash
bash /data/pipelines/ncov2019-artic-nf/tests/GPAS_tests.bash
```

This should run through each sequencing technology and bioinformatics workflow combination. It will report any discrepencies between the outputs produced and outputs expected.

## Test dataset

For this testing the "core" test samples are used. These are

### Illumina
```bash
SARS-CoV-2_reference_ox,mmm-artic-ill-s11511-1
SARS-CoV-2_reference_ox,mmm-artic-ill-s12220-1
SARS-CoV-2_reference_ox,mmm-artic-ill-s12368-1
SARS-CoV-2_reference_ox,mmm-artic-ill-s16621-1
SARS-CoV-2_reference_ox,mmm-artic-ill-s24350-3
SARS-CoV-2_reference_ox,mmm-artic-ill-s32219-2
SARS-CoV-2_reference_ox,mmm-artic-ill-s53667-1
SARS-CoV-2_reference_ox,mmm-artic-ill-s59130-3
SARS-CoV-2_reference_ox,mmm-artic-ill-s64379-3
SARS-CoV-2_reference_ox,mmm-artic-ill-s71898-2
SARS-CoV-2_reference_ox,mmm-artic-ill-s82718-1
SARS-CoV-2_reference_ox,mmm-artic-ill-s98244-3
```

### ONT
```bash
SARS-CoV-2_reference_ox,mmm-artic-ont-s11511-1
SARS-CoV-2_reference_ox,mmm-artic-ont-s12220-4
SARS-CoV-2_reference_ox,mmm-artic-ont-s12368-1
SARS-CoV-2_reference_ox,mmm-artic-ont-s16621-3
SARS-CoV-2_reference_ox,mmm-artic-ont-s24350-1
SARS-CoV-2_reference_ox,mmm-artic-ont-s32219-1
SARS-CoV-2_reference_ox,mmm-artic-ont-s53667-1
SARS-CoV-2_reference_ox,mmm-artic-ont-s59130-1
SARS-CoV-2_reference_ox,mmm-artic-ont-s64379-1
SARS-CoV-2_reference_ox,mmm-artic-ont-s71898-1
SARS-CoV-2_reference_ox,mmm-artic-ont-s82718-2
SARS-CoV-2_reference_ox,mmm-artic-ont-s98244-1
```

## Outputs

The outputs for each run are stored in `/work/outputs/`. There will be four folders for these tests:

* /work/output/illumina_viridian_test/
* /work/output/illumina_artic_test/
* /work/output/ont_viridian_test/
* /work/output/ont_artic_test/

Within these folders, the tests summarise the outputs into a tsv file, for example `illumina_viridian_summary.tsv`. There is also a comparison file containing any discrepencies in `illumina_viridian_comparison.tsv`.
