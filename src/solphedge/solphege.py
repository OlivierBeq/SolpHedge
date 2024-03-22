# -*- coding: utf-8

"""pH-dependent solubility machine learning model."""

import json
import os
import pickle
import sys
import warnings

import ml2json
import numpy as np
import pandas as pd
from BlueDesc_pywrapper import BlueDesc
from CDK_pywrapper import CDK
from Mold2_pywrapper import Mold2
from PaDEL_pywrapper import PaDEL, descriptors as PaDEL_descriptors
from chembl_structure_pipeline import standardize_mol as csp_standardize
from chemopy import ChemoPy
from papyrus_structure_pipeline import standardize as psp_standardize
from rdkit import Chem


# Filter out warnings of ml2json about scikit-learn's version when loading the scalers
warnings.filterwarnings(action='ignore', category=UserWarning)
warnings.filterwarnings(action='ignore', category=FutureWarning)


class NumPyMLPRegressor:
    '''
        Multilayer Perceptron for regression.
        '''

    def __init__(self, weight_file: str):
        """Load a fully connected neural network model for inference.

        :param weight_file: file on disk storing the model's weights
        """
        # Read pickle file
        with open(weight_file, 'rb') as fh:
            self.weights = pickle.load(fh)

    def predict(self, x: np.ndarray):
        """Run inference on the input

        :param x: Input data
        :return: Predictions
        """
        for i in range(1, 5):
            # Weights
            x = np.dot(x, self.weights[f'fc{i}.weight'].T)
            # Bias
            x += self.weights[f'fc{i}.bias']
            # ReLU
            x = (abs(x) + x) / 2
        return x


class SolpH:
    """Machine learning model for pH-dependent predictions."""

    path_fn = (lambda x: os.path.abspath(os.path.join(os.path.dirname(__file__), 'models', x)))

    # Default models
    _modelfiles = {
        # Model predicting solubility (negative log scale) at pH=7.4
        'negLogS_ph7.4': {'scaler': path_fn('neglogS_scaler.json'),
                           'novariance': path_fn('neglogS_features_novariance.json'),
                           'model': path_fn('neglogS_FFNN_model.pkl')
                           },
        # Model predicting the difference in solubility between pH=7.4 and pH=1.0
        'deltaLogS': {'scaler': path_fn('deltalogS_scaler.json'),
                       'novariance': path_fn('deltalogS_features_novariance.json'),
                       'model': path_fn('deltalogS_FFNN_model.pkl')
                       }
    }

    def __init__(self, standardize: bool = True, standardizer: str = 'papyrus', njobs: int = 1):
        """Instantiate the machine learning model for prediction of pH-dependent solubility.

        :param standardize: should molecules be standardized prior to predictions
        :param standardizer: standardizer to apply; one of {papyrus, chembl}
        """
        if sys.platform not in ['win32', 'linux']:
            raise RuntimeError(f'SolpH can only be used on Windows and Linux platforms.')
        if standardizer.lower() not in ['papyrus', 'chembl']:
            raise ValueError('Standardizer must be either \'papyrus\' or \'chembl\'.')
        self.standardize = standardize
        self.standardizer = standardizer.lower()
        self.njobs = njobs
        # Instantiate models
        #       1) deltaLogS
        self._deltaLogS_scaler = ml2json.from_json(self._modelfiles['deltaLogS']['scaler'])
        with open(self._modelfiles['deltaLogS']['novariance']) as fh:
            self._deltaLogS_novar_features = json.load(fh)
        self._deltaLogS_model = NumPyMLPRegressor(self._modelfiles['deltaLogS']['model'])
        deltaLogS_n_inputs = self._deltaLogS_model.weights['fc1.weight'].shape[1]
        if self._deltaLogS_scaler.n_features_in_ - len(self._deltaLogS_novar_features) != deltaLogS_n_inputs:
            raise ValueError('Model files for deltaLogS do not match one another. Contact the maintainer.')
        #       2) negLogS
        self._negLogS_scaler = ml2json.from_json(self._modelfiles['negLogS_ph7.4']['scaler'])
        with open(self._modelfiles['negLogS_ph7.4']['novariance']) as fh:
            self._negLogS_novar_features = json.load(fh)
        self._negLogS_model = NumPyMLPRegressor(self._modelfiles['negLogS_ph7.4']['model'])
        negLogS_n_inputs = self._negLogS_model.weights['fc1.weight'].shape[1]
        if self._negLogS_scaler.n_features_in_ - len(self._negLogS_novar_features) != negLogS_n_inputs:
            raise ValueError('Model files for negLogS do not match one another. Contact the maintainer.')


    def predict(self, mols: Chem.Mol | list[Chem.Mol] | str | list[str], round: int = 3, out_units: str='mM'):
        """Predict delta solubility values between pH of 7.4 and 1.0.

        :param mols: RDKit molecule(s) or SMILES
        :param round: number of decimals to keep in the predictions
        :param out_units: units of molar concentrations {M; mM; uM; nM; pM}
        """
        allowed_units = {'M': 0, 'mM': -3, 'uM': -6, 'nM': -9, 'pM': -12}
        if out_units not in allowed_units.keys():
            raise ValueError(f'Output units must be one of {allowed_units}')
        smiles = None
        # Ensure a list is processed if only one input provided
        if isinstance(mols, (Chem.Mol, str)):
            mols = [mols]
        # If SMILES provided, parse them
        if isinstance(mols, str) or all(isinstance(mol, str) for mol in mols):
            smiles = mols[:]
            mols = [Chem.MolFromSmiles(x) for x in smiles]
        # Raise error neither SMILES not rdkit molecules are provided
        if not isinstance(mols, (Chem.Mol, type(None))) and not all(isinstance(mol, (Chem.Mol, type(None))) for mol in mols):
            raise ValueError('Can only predict solubilities for RDKit molecule object(s) or SMILES string(s).')
        # If only one molecule provided
        # If rdkit molecules provided, obtain SMILES
        if smiles is None:
            smiles = [Chem.MolToSmiles(x) for x in mols]
        # Standardize if need be:
        if self.standardize:
            if self.standardizer == 'papyrus':
                mols = [psp_standardize(mol, raise_error=False,
                                        filter_non_small_molecule=False,
                                        filter_mixtures=False,
                                        filter_inorganic=False) for mol in mols]
            elif self.standardizer == 'chembl':
                mols = [csp_standardize(mol) for mol in mols]
            else:
                raise NotImplementedError(f'Standardizer \'{self.standardizer}\' not implemented.')
        # Identify None values
        none_idx = np.where(np.equal(mols, None))[0]
        mols = [mol for mol in mols if mol is not None]
        # Add hydrogens
        mols = [Chem.AddHs(mol) for mol in mols]
        # Calculate Mold2 descriptors
        mold2 = Mold2(fill_na=0, verbose=False)
        descs_mold2 = mold2.calculate(mols, show_banner=False, njobs=self.njobs, chunksize=-(-len(mols) // self.njobs))
        descs_mold2.columns = 'Mold2_' + descs_mold2.columns
        chemopy = ChemoPy()
        descs_chemopy = chemopy.calculate(mols, show_banner=False, njobs=self.njobs, chunksize=-(-len(mols) // self.njobs))
        descs_chemopy.columns = 'ChemoPy_' + descs_chemopy.columns
        padel = PaDEL(PaDEL_descriptors)
        descs_padel = padel.calculate(mols, show_banner=False, njobs=self.njobs, chunksize=-(-len(mols) // self.njobs))
        descs_padel.columns = 'PaDEL_' + descs_padel.columns
        cdk = CDK()
        descs_cdk = cdk.calculate(mols, show_banner=False, njobs=self.njobs, chunksize=-(-len(mols) // self.njobs))
        descs_cdk.columns = 'CDK_' + descs_cdk.columns
        bluedesc = BlueDesc()
        descs_bluedesc = bluedesc.calculate(mols, show_banner=False, njobs=self.njobs, chunksize=-(-len(mols) // self.njobs))
        descs_bluedesc.columns = 'BlueDesc_' + descs_bluedesc.columns
        # Combine datasets
        data = pd.concat([descs_mold2, descs_chemopy, descs_padel, descs_cdk, descs_bluedesc], axis=1)
        del descs_mold2, descs_chemopy, descs_padel, descs_cdk, descs_bluedesc
        # Scale the data
        negLogS_data = pd.DataFrame(self._negLogS_scaler.transform(data), columns=data.columns)
        deltaLogS_data = pd.DataFrame(self._deltaLogS_scaler.transform(data), columns=data.columns)
        del data
        # Drop columns with low variance
        negLogS_data = negLogS_data.drop(columns=self._negLogS_novar_features)
        deltaLogS_data = deltaLogS_data.drop(columns=self._deltaLogS_novar_features)
        # Predict the solubility difference
        negSol = self._negLogS_model.predict(negLogS_data.values)
        deltaSol = self._deltaLogS_model.predict(deltaLogS_data.values)
        # Insert missing values where needed
        negSol = np.insert(negSol, none_idx, np.NaN, axis=0).ravel()
        deltaSol = np.insert(deltaSol, none_idx, np.NaN, axis=0).ravel()
        # Combine predictions
        preds = pd.concat([pd.Series(smiles, name='molecule'),
                           pd.Series(negSol, name='-logS (pH=7.4)'),
                           pd.Series(deltaSol, name='delta logS (pH=7.4 - pH=1.0)'),
                           pd.Series(negSol + deltaSol, name='composite -logS (pH=1.0)')],
                          axis=1)
        preds[f'solubility (pH=7.4; {out_units})'] = 10 ** -(preds['-logS (pH=7.4)'] + allowed_units[out_units])
        preds[f'solubility (pH=1.0; {out_units})'] = 10 ** -(preds['composite -logS (pH=1.0)'] + allowed_units[out_units])
        # Round the data
        return preds.round(round)
