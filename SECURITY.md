# Safety, scope, and reporting

## Scope of this project

This is exploratory computational research software. It is **not** medical
advice, **not** a diagnostic device, **not** a treatment recommendation, and
**not** evidence of efficacy or safety for any multiple-sclerosis intervention.

Do not use any output from this repository to make decisions about a person's
care, to enrol or exclude trial participants, to change medication, or as a
substitute for qualified clinicians, clinical trials, regulatory review, or
published evidence.

See [`docs/QUALITY.md`](docs/QUALITY.md) for what the repo's own labels claim
and, more importantly, what they do not.

## Data

No patient data is in this repository and none should ever be added. The
loaders fetch open published datasets at runtime into `data/cache/`, which is
git-ignored. No model weights and no datasets are committed or distributed here.

## Reporting a vulnerability

Report security or privacy issues privately to the repository owner via GitHub
Security Advisories on this repository, rather than in a public issue.

Never include patient data, private data, or credentials in a report, an issue,
or a pull request. A minimal reproducible description is enough.
