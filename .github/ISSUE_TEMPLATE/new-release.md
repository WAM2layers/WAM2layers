---
name: New release
about: Release checklist
title: Release vX.Y.Z
labels: ''
assignees: '' 
---

## Release checklist
- [ ] Create a new branch, e.g. `prep_release_vXXX`
- [ ] Update "version" in `src/wam2layers/__init__.py` (use [semantic versioning](https://semver.org/))
- [ ] Make sure the changelog is up to date; update the "Unreleased" header with the
  new version, and add a new "Unreleased" header above
- [ ] Make sure the code is neatly formatted according to the [style guide](https://wam2layers.readthedocs.io/en/latest/develop.html#follow-the-style-guide).
- [ ] Review and merge the release branch
- [ ] Check on Zenodo whether the Zenodo integration with WAM2layers/WAM2layers is switched to ON. In case the switch is OFF, do not turn it ON, but contact @ruudvdent to turn it ON for you. This is in order to avoid creating multiple Webhooks and thus multiple releases.
- [ ] [Make a release on GitHub](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)
- [ ] Upon release, the [Github Actions workflow for publishing](https://github.com/WAM2layers/WAM2layers/actions/workflows/publish.yaml) will run. More extensive documentation on publishing is available [here](https://packaging.python.org/en/latest/guides/publishing-package-distribution-releases-using-github-actions-ci-cd-workflows/).
- [ ] Check that the new release is successfully rendered on Zenodo and PyPI.
- [ ] Verify that you can install the new version with pip.
- [ ] Celebrate.


---
name: New release
about: Release checklist
title: Release 2025.X.Y
labels: ''
assignees: '' 
---
