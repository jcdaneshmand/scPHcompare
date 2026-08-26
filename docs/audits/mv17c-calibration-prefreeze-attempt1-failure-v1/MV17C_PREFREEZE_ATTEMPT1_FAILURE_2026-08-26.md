# MV17-C prefreeze attempt 1 classification stop

Attempt 1 at implementation head `5725915` stopped before writing any private
locator or public audit artifact and before selecting or generating a null.
The builder incorrectly required the MV08-Q newly produced source table to have
132 rows; that closure correctly contains 129 new sources. The accepted MV08-R
binding composes those 129 sources with three previously qualified sources for
the required 132-unit axis.

Recovery changes only the source-binding classification: use the exact 132-row
MV08-R binding, add its manifest to the prerequisite rehash, and admit the fifth
prior-source cache root. Sentinel selection, null families, seeds, summaries,
thresholds, resources, privacy rules, and downstream firewalls are unchanged.
Attempt-1 roots are preserved and must not be reused.
