# Source baseline

`chip2tracks` began as a clean working-tree snapshot of local
`cutnrun2tracks` version 0.2.8 on 2026-08-25. The source repository had no
commits and no configured remote, so a commit hash could not provide a
trustworthy identity.

The file-level SHA-256 manifest is
`provenance/cutnrun2tracks_v0.2.8.sha256`. It contains 84 files, excludes only
the nested `.git` directory, and has SHA-256 digest:

`4ccc8c329a52a33803e1425149b17eeb55751c74cfa2410704c210f89e5b8121`

Recommended history before publication:

1. Commit the unmodified import and tag it
   `cutnrun2tracks-source-v0.2.8-local`.
2. Apply the `chip2tracks` changes as later commits.
3. Never rewrite or move the baseline tag.
4. Preserve this manifest in releases so the local source can be verified
   independently of Git history.

This is stronger than claiming equivalence to a GitHub tag that was not
actually present in the local repository.

Verify the still-available local source snapshot with:

```bash
bash utilities/verify_cutnrun_baseline.sh ../cutnrun2tracks
```
