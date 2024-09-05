

# <sup>:notes:</sup>:saxophone: Sol*pH*<sup>edge</sup>

pH-dependent solubility predictions for small molecules

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## This branch describes the traning process of models

:warning: MacOS is not supported by SolpH<sup>edge</sup>. :warning:

![Windows](https://img.shields.io/badge/Windows-0078D6?style=for-the-badge&logo=windows&logoColor=white) ![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)

### 1) Dependencies
The easiest way to train SolpH<sup>edge</sup> models is to create a new conda environment:

```bash
conda env create --file=requirements.yaml
```

You may also inspect the file and install dependencies yourself. 

### 2) Training

The python script called `training.py` trains SolpH<sup>edge</sup> models and can be used like so:

```bash
cd example
python ../solphedge_training.py train -i LogS.tsv -e LogS -o LogS_FCNN --split random --smilesColumn SMILES 
```

### 3) Exporting PyTorch models for NumPy

To ensure minimal dependency on python libraries at inference, models' state dictionaries are exported to pickle files of dictionaries of NumPy arrays.

This step is carried out with:

```bash
python solphedge_training.py convert -i LogS_FCNN.pkg -o LogS_FCNN_numpy
```
