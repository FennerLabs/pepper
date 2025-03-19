# PEPPER - Predict Environmental Pollutant PERsistence

PEPPER is a package developed by the Fenner Labs for analyzing and modeling persistence of micropollutants in different environments.

## Installation 

PEPPER may be installed using:
```
pip install pepper
```

or clone this repository
```
git clone https://github.com/FennerLabs/pepper
git cd pepper
```

We also recommend creating a dedicated virtual environment with python 3.9 as base 

```
python -m venv pepper_env
source pepper_env/bin/activate
```

We have included all requirements in the pyproject.toml file so all dependencies may be installed as follows

``` 
pip install .
```

After installation you may use all pepper methods and reproduceworkflows from our current projects

## Projects

### Bayesian inference for soil biotransformation half-lives - Hafner et al., 2023
Here's how to reproduce the data and the figures from the publication:
```
cd scripts
python bayesian_inference_main.py
```

### Current Opportunities and Limitations in Predicting Micropollutant Removal in Wastewater Treatment based on Molecular Structure - Cordero et al., 2025
In this project we include methods to model the breakthrough of micropollutants in wastewater treatment plants.
Main results can be reproduced as follows:
```
cd scripts
predict_breakthrough_wwtp.py
```

Please refer to the main publication for further details

## Sessions 
Use this link to start a session and test PEPPER
[![launch - renku](https://renkulab.io/renku-badge.svg)](https://renkulab.io/projects/fenner-labs/projects/pepper/sessions/new?autostart=1)

## Related Projects
We also have a pepper_app 
[<img alt="launch - streamlit" height="20" src="https://streamlit.io/images/brand/streamlit-mark-color.svg" title="Launch pepper_app" width="20"/>](https://pepper-app.streamlit.app)
to predict several endpoints of interest related to environmental persistence