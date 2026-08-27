# Test matrix

The synthetic suite covers Bash syntax, safe config parsing, both assay enums,
optional/exact/shared controls, narrow/broad caller routing, control-free SE,
peak-caller continuation, filtering cleanup scope, consensus zero-count
tolerance, genome-order annotation, reporting, source-baseline format and
scientific helper functions.

CI runs:

```bash
python -m compileall -q scripts common tests
bash tests/run_tests.sh
```

Before release, add small real datasets for:

- PE narrow ChIP-seq with two replicates and matched inputs;
- broad ChIP-seq with epic2 and MACS3 broad comparison;
- PE ChIPmentation with adapter auto-detection and an explicit override;
- SE ChIP-seq fragment modeling;
- control-free peak calling with explicit report warning;
- competitive host+dm6 spike-in, including low-count warning behavior;
- default two-biological-sample consensus, optional pairwise IDR, broad
  consensus, FRiP, NRF/PBC, NSC/RSC, PCA and correlation;
- restart from each scientifically meaningful stage.

Freeze expected ranges rather than exact bytes for stochastic or
version-sensitive biological outputs. A release candidate fails if its values
leave the justified ranges or if an expected warning disappears.
