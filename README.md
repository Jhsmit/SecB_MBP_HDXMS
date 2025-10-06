# How SecB maintains clients in a translocation competent state

This repository contains scripts to reproduce HDX-MS results from the paper "[How SecB maintains clients in a translocation competent state](https://doi.org/10.1038/s42003-025-08821-2)"


Datasets are published in hdxms-datasets v0.1.5 (legacy) format in the `data` directory. 


## Install


Create an environment:

```
uv venv
```
Activate the environment, then install the requirements:

```
uv pip install -r requirements.txt
```

Install PyHDX / smitfit from GitHub:

```
uv pip install git+https://github.com/Jhsmit/pyhdx@ce183307df114812611c9ad8e7bf7053e7dd43bd"
```

```
uv pip install git+https://github.com/Jhsmit/smitfit@891cb68f229d372e352b7612106a65a5f7129140"
```


Generate requirements from `freeze.txt` (optional):

```
findstr /v "^-e " freeze.txt > requirements.txt
```
or
```
grep -v "^-e " freeze.txt > requirements.txt
```

