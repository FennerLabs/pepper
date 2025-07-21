# Changelog

## [Unreleased]
- FIX: Incomplete config files and scripts which did not include all regressors were accidentally uploaded. I have added the new versions and confirmed that these work regardless of the operating system. 
- update: changed mean_squared_error(squared=False) to use root_mean_squared_error as this is the new recommended way. 