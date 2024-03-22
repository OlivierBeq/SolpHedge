[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

![Windows](https://img.shields.io/badge/Windows-0078D6?style=for-the-badge&logo=windows&logoColor=white) ![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)  

# <sup>:notes:</sup>:saxophone: Sol*pH*<sup>edge</sup>

pH-dependent solubility predictions for small molecules


## Installation

:warning: MacOS is not supported by SolpH<sup>edge</sup>. :warning:

OpenBabel is the only dependency of SolpH<sup>edge</sup> not available through PyPI. It can be installed with:

```bash
conda install openbabel -c conda-forge
```

Then SolpH<sup>edge</sup> and other dependencies can be isntalled with: 

```bash
pip install solphedge
```

## Get started

```python
from solphedge import SolpH

smiles_list = ['CC(C1=CC(=CC=C1)C(=O)C2=CC=CC=C2)C(=O)O',
               'CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O',
               'C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl',
               'C1=CC(=C(C=C1C2=C(C=C(C=C2)F)F)C(=O)O)O',
               'CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC',
               'C1=CC=C(C=C1)CC2NC3=C(C=C(C(=C3)C(F)(F)F)S(=O)(=O)N)S(=O)(=O)N2',
               'CC1=CN=C(S1)NC(=O)C2=C(C3=CC=CC=C3S(=O)(=O)N2C)O',
               'C#CCO[C@H]1CN2CCC1CC2',
               'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O']
model = SolpH()
print(model.predict(smiles_list))
#                                             molecule  -logS (pH=7.4)  delta logS (pH=7.4 - pH=1.0)  composite -logS (pH=1.0)  solubility (pH=7.4; mM)  solubility (pH=1.0; mM)
# 0            CC(C1=CC(=CC=C1)C(=O)C2=CC=CC=C2)C(=O)O           3.859                         2.910                     6.770                    0.138                    0.000
# 1       CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O           4.793                         0.578                     5.370                    0.016                    0.004
# 2           C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl           3.792                         3.072                     6.864                    0.161                    0.000
# 3            C1=CC(=C(C=C1C2=C(C=C(C=C2)F)F)C(=O)O)O           4.113                         2.768                     6.881                    0.077                    0.000
# 4   CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC           3.862                         0.000                     3.862                    0.137                    0.137
# 5  C1=CC=C(C=C1)CC2NC3=C(C=C(C(=C3)C(F)(F)F)S(=O)...           4.042                         0.127                     4.169                    0.091                    0.068
# 6   CC1=CN=C(S1)NC(=O)C2=C(C3=CC=CC=C3S(=O)(=O)N2C)O           4.723                         0.000                     4.723                    0.019                    0.019
# 7                              C#CCO[C@H]1CN2CCC1CC2           3.749                         0.000                     3.749                    0.178                    0.178
# 8                      CC(C)CC1=CC=C(C=C1)C(C)C(=O)O           3.920                         1.473                     5.393                    0.120                    0.004
```

## Other parameters

```python
class SolpH(standardize: bool = True, standardizer: str = 'papyrus', njobs: int = 1):
    ...
```

#### Parameters

- ***standardize  : bool***  
  Should input molecules or SMILES be standardized before making predictions.
- ***standardizer  : str = {['chembl'](https://github.com/chembl/ChEMBL_Structure_Pipeline); ['papyrus'](https://github.com/OlivierBeq/Papyrus_structure_pipeline)}***  
  Molecular standardizer to be applied.
- ***njobs  : int***  
  Maximum number of simultaneous processes calculating molecular descriptors prior to predictions.

```python
    ...
    def predict(mols: Chem.Mol | list[Chem.Mol] | str | list[str], round: int = 3, out_units: str='mM'):
```

#### Parameters

- ***mols      : str | rdkit.Chem.Mol | list[str] | list[rdkit.Chem.Mol]***  
  Input molecule(s) or SMILES.
- ***round     : int***  
  Number of decimals to round the result.
- ***out_units : str = {'M'; 'mM'; 'uM'; 'nM'; 'pM'}***  
  Units to convert log-scale molar solubilities to.
