# -*- coding: utf-8 -*-

import os
import json
import pickle
import time
from collections import defaultdict
from datetime import datetime, timedelta

import click
import numpy as np
import pandas as pd
import ml2json

from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, MaxAbsScaler
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import pearsonr
from tqdm import tqdm

from Mold2_pywrapper import Mold2
from chemopy import ChemoPy
from PaDEL_pywrapper import PaDEL, descriptors as PaDEL_descriptors
from CDK_pywrapper import CDK
from BlueDesc_pywrapper import BlueDesc

import torch
from torch import nn

import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

PATIENCE = 500
MIN_EPOCHS = 1000
LEARNING_RATE = 1e-4
SEED = 12345678

raw_input = input


def generate_scaffold(smiles, include_chirality=False):
    """return scaffold string of target molecule"""
    mol = Chem.MolFromSmiles(smiles)
    scaffold = MurckoScaffoldSmiles(mol=mol, includeChirality=include_chirality)
    return scaffold


class ScaffoldSplitter:
    """Class for doing data splits by chemical scaffold.

    Referr to Deepchem for the implementation, https://git.io/fXzF4
    """

    def _split(self, dataset: pd.DataFrame,
               smiles_list: list[str],
               frac_train: float = 0.8,
               frac_valid: float = 0.1,
               frac_test: float = 0.1,
               random_state: int = 0,
               include_chirality: bool = False):
        if not isinstance(dataset, pd.DataFrame):
            raise ValueError("dataset must be a pandas dataframe")
        np.testing.assert_almost_equal(frac_train + frac_valid + frac_test, 1.)
        if len(dataset) != len(smiles_list):
            raise ValueError(
                f"The lengths of dataset({len(dataset)}) and smiles_list ({len(smiles_list)}) are different")
        smiles_list = list(smiles_list)

        rng = np.random.RandomState(random_state)

        scaffolds = defaultdict(list)
        for ind, smiles in enumerate(tqdm(smiles_list, desc='Obtaining scaffolds')):
            scaffold = generate_scaffold(smiles, include_chirality)
            scaffolds[scaffold].append(ind)

        scaffold_inds = list(scaffolds.values())
        scaffold_sets = [scaffold_inds[i] for i in rng.permutation(list(range(len(scaffold_inds))))]

        n_total_valid = int(np.floor(frac_valid * len(dataset)))
        n_total_test = int(np.floor(frac_test * len(dataset)))

        train_index = []
        valid_index = []
        test_index = []

        for scaffold_set in tqdm(scaffold_sets, desc='Gathering scaffolds'):
            if len(valid_index) + len(scaffold_set) <= n_total_valid:
                valid_index.extend(scaffold_set)
            elif len(test_index) + len(scaffold_set) <= n_total_test:
                test_index.extend(scaffold_set)
            else:
                train_index.extend(scaffold_set)

        return np.array(train_index), np.array(valid_index), np.array(test_index)

    def train_valid_test_split(self, dataset: pd.DataFrame,
                               smiles_list: list[str],
                               frac_train: float = 0.8,
                               frac_valid: float = 0.1,
                               frac_test: float = 0.1,
                               random_state: int = 0,
                               include_chirality: bool = False,
                               verbose: bool = True,
                               logfile=None):
        """Split dataset into train, valid and test set.

        Split indices are generated by splitting based on the scaffold of small
        molecules.

        :param dataset: pandas dataframe
        :param smiles_list: SMILES list corresponding to the dataset
        :param frac_train: fraction of dataset put into training data
        :param frac_valid: fraction of dataset put into validation data
        :param frac_test: fraction of dataset put into test data
        :param random_state: random seed
        :param include_chirality: should the scaffolds include chirality
        :param verbose: display splitting information
        :return: indices of the dataset split into training, validation and test
        """
        train_inds, valid_inds, test_inds = self._split(dataset,
                                                        smiles_list,
                                                        frac_train,
                                                        frac_valid,
                                                        frac_test,
                                                        random_state,
                                                        include_chirality)

        if verbose:
            self.log_results(len(dataset), (frac_train, frac_valid, frac_test),
                             len(train_inds), len(valid_inds), len(test_inds),
                             logfile=logfile)

        return train_inds, valid_inds, test_inds

    def train_valid_split(self, dataset: pd.DataFrame,
                          smiles_list: list[int],
                          frac_train: float = 0.9,
                          frac_valid: float = 0.1,
                          random_state: int = 0,
                          include_chirality: bool = False,
                          verbose: bool = True,
                          logfile=None):
        """Split dataset into train and valid set.

        Split indices are generated by splitting based on the scaffold of small
        molecules.

        :param dataset: pandas dataframe
        :param smiles_list: SMILES list corresponding to the dataset
        :param frac_train: fraction of dataset put into training data
        :param frac_valid: fraction of dataset put into validation data
        :param random_state: random seed
        :param include_chirality: should the scaffolds include chirality
        :param verbose: display splitting information
        :return: indices of the dataset split into training and validation
        """
        train_inds, valid_inds, test_inds = self._split(dataset, smiles_list,
                                                        frac_train,
                                                        frac_valid,
                                                        0.,
                                                        random_state,
                                                        include_chirality)
        assert len(test_inds) == 0

        if verbose:
            self.log_results(len(dataset), (frac_train, frac_valid),
                             len(train_inds), len(valid_inds), logfile=logfile)

        return train_inds, valid_inds

    def log_results(self, total_len: int,
                    ideal_ratios: tuple[float, float] | tuple[float, float, float],
                    len_train_inds: int, len_valid_inds: int, len_test_inds: int = None,
                    logfile=None):
        ideal_train_len = int(round(ideal_ratios[0] * total_len))
        train_error = (len_train_inds - ideal_train_len) / ideal_train_len
        print(f'#rows in taining set: {len_train_inds}, ideally: {ideal_train_len}, error: {train_error:.5%}',
              file=logfile, flush=True)
        ideal_valid_len = int(round(ideal_ratios[1] * total_len))
        valid_error = (len_valid_inds - ideal_valid_len) / ideal_valid_len
        print(f'#rows in validation set: {len_valid_inds}, ideally: {ideal_valid_len}, error: {valid_error:.5%}',
              file=logfile, flush=True)
        if len_test_inds is not None:
            ideal_test_len = int(round(ideal_ratios[2] * total_len))
            test_error = (len_test_inds - ideal_test_len) / ideal_test_len
            print(f'#rows in validation set: {len_test_inds}, ideally: {ideal_test_len}, error: {test_error:.5%}',
                  file=logfile, flush=True)


class MLPRegressor(nn.Module):
    '''
    Multilayer Perceptron for regression.
    '''

    def __init__(self, in_features: int):
        super().__init__()
        self.fc1 = nn.Linear(in_features, 1000)
        self.fc2 = nn.Linear(1000, 200)
        self.fc3 = nn.Linear(200, 100)
        self.fc4 = nn.Linear(100, 1)
        self.act = nn.functional.relu
        self.dropout = nn.Dropout(0.1)

    def forward(self, x, dropout: bool = False):
        '''
            Forward pass
        '''
        res = self.act(self.fc1(x))
        if dropout:
            res = self.dropout(res)
        res = self.act(self.fc2(res))
        if dropout:
            res = self.dropout(res)
        res = self.act(self.fc3(res))
        if dropout:
            res = self.dropout(res)
        res = self.fc4(res)
        return res

    @property
    def num_params(self):
        model_parameters = filter(lambda p: p.requires_grad, self.parameters())
        params = sum([np.prod(p.size()) for p in model_parameters])
        return params


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    """Group allowing subcommands to be defined"""
    pass

@main.command(help='Train a SolpHedge model.', context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input', type=str, required=True, multiple=False,
              default=None, nargs=1, show_default=False,
              help='Input tab-separated file with extension containing the molecular SMILES and '
                   'solubility values to train the model from.')
@click.option('-e', '--endpoint', type=click.Choice(['LogS', 'deltaLogS'], case_sensitive=False),
              required=True, default=None, nargs=1,
              show_default=False, multiple=False,
              help=('Endpoint to train the model from. Must be either \'LogS\' or \'deltaLogS\'. '
                    'The variable must appear as a data column in the input file.'))
@click.option('-o', '--output', type=str, required=True, multiple=False,
              default=None, nargs=1, show_default=False,
              help='Common name of the output files without extension '
                   '(weights of PyTorch model, feature scaler, dropped feature names, logfile).')
@click.option('--split', type=click.Choice(['random', 'scaffold'], case_sensitive=False),
              required=False, default='random', nargs=1, show_default=True,
              help='Type of data splitting scheme to obtain training, validation and test sets.')
@click.option('--smilesColumn', type=str, required=False, default='SMILES', nargs=1,
              show_default=True, help='Name of the data column containing SMILES values.')
def train(input: str, endpoint: str, output: str, split: str, smilescolumn: str):
    """CLI to train SolpHedge models."""
    endpoint = endpoint.lower()
    split = split.lower()
    smilescolumn = smilescolumn.lower()
    if endpoint not in ['logs', 'deltalogs']:
        raise ValueError('Endpoint must be either \'LogS\' or \'deltaLogS\'.')
    # Ensure input file exists
    if not os.path.exists(input):
        raise FileNotFoundError(f'File does not exist: {input}')
    # Ensure output file does not exist or ask if overwriting
    if (os.path.exists(f'{output}_scaler.json') or
            os.path.exists(f'{output}_dropped_features.json') or
            os.path.exists(f'{output}_logfile.log') or
            os.path.exists(f'{output}_model.pkg')):
        while (confirm := raw_input('Output file(s) already exists. Overwrite? (y/N)\n')).lower() not in ['y', 'n']:
            pass
        if confirm.lower() == 'n':
            return
    # Ensure columns exist in the input file
    snippet = pd.read_csv(input, sep='\t', nrows=100)
    snippet.columns = snippet.columns.str.lower()
    if endpoint not in snippet.columns:
        raise ValueError(f'Missing column {endpoint} in input file.')
    if smilescolumn not in snippet.columns:
        raise ValueError(f'Missing column {smilescolumn} in input file.')
    del snippet
    # Create log file
    logfile = open(f'{output}_logfile.log', 'w')
    start_time = time.monotonic()
    logfile.write(f'Start time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
    # Read entire file
    data = pd.read_csv(input, sep='\t')
    data.columns = data.columns.str.lower()
    # Parse molecules
    data['mol'] = [Chem.MolFromSmiles(smiles) for smiles in data[smilescolumn]]
    # Log the number of unparseable molecules
    if data['mol'].isnull().any():
        print(f'{data["mol"].isnull().sum()} molecules could not be parsed.\n')
        data = data.dropna(subset='mol')
    # Add hydrogen atoms
    data['mol'] = data['mol'].apply(Chem.AddHs)
    # Drop duplicate molecules based on connectivity
    data['connectivity'] = data['mol'].apply(lambda x: Chem.MolToInchiKey(x).split('-')[0])
    n_mols_orig = len(data)
    data = data.drop_duplicates(subset='connectivity').reset_index(drop=True)
    if n_mols_orig != len(data):
        print(f'{n_mols_orig - len(data)} duplicated molecules (2D topology) have been dropped.\n'
              f'The model will be trained from the remaining {len(data)}.')
    # Calculate 2D molecular descriptors
    mold2_descs = Mold2(verbose=False).calculate(data['mol'], show_banner=False, njobs=-1, chunksize=100)
    # Rename and deal with abnormal values
    mold2_descs.rename(columns={colname: f'Mold2_{colname}' for colname in mold2_descs.columns}, inplace=True)
    mold2_descs = mold2_descs.fillna(0)
    mold2_descs = (mold2_descs
                   .apply(pd.to_numeric, downcast='integer')
                   .apply(pd.to_numeric, downcast='float'))
    mold2_descs[mold2_descs.abs() > 1e5] = 0
    chemo_descs = ChemoPy().calculate(data['mol'], show_banner=False, njobs=-1, chunksize=100)
    chemo_descs.rename(columns={colname: f'ChemoPy_{colname}' for colname in chemo_descs.columns}, inplace=True)
    chemo_descs = chemo_descs.fillna(0)
    chemo_descs = (chemo_descs
                   .apply(pd.to_numeric, downcast='integer')
                   .apply(pd.to_numeric, downcast='float'))
    chemo_descs[chemo_descs.abs() > 1e5] = 0
    padel_descs = PaDEL(PaDEL_descriptors).calculate(data['mol'], show_banner=False, njobs=-1, chunksize=100)
    padel_descs.rename(columns={colname: f'PaDEL_{colname}' for colname in padel_descs.columns}, inplace=True)
    padel_descs = padel_descs.fillna(0)
    padel_descs[padel_descs.abs() > 1e5] = 0
    padel_descs = (padel_descs
                   .apply(pd.to_numeric, downcast='integer')
                   .apply(pd.to_numeric, downcast='float'))
    cdk_descs   = CDK().calculate(data['mol'], show_banner=False, njobs=-1, chunksize=100)
    cdk_descs.rename(columns={colname: f'CDK_{colname}' for colname in cdk_descs.columns}, inplace=True)
    cdk_descs = cdk_descs.fillna(0)
    cdk_descs[cdk_descs.abs() > 1e5] = 0
    cdk_descs = (cdk_descs
                 .apply(pd.to_numeric, downcast='integer')
                 .apply(pd.to_numeric, downcast='float'))
    blued_descs = BlueDesc().calculate(data['mol'], show_banner=False, njobs=-1, chunksize=100)
    blued_descs.rename(columns={colname: f'BlueDesc_{colname}' for colname in blued_descs.columns}, inplace=True)
    blued_descs = blued_descs.fillna(0)
    blued_descs[blued_descs.abs() > 1e5] = 0
    blued_descs = (blued_descs
                   .apply(pd.to_numeric, downcast='integer')
                   .apply(pd.to_numeric, downcast='float'))
    blued_descs[blued_descs.abs() > 1e5] = 0
    # Combine descriptors
    descs = pd.concat((mold2_descs, chemo_descs, padel_descs, cdk_descs, blued_descs), axis=1)
    del mold2_descs, chemo_descs, padel_descs, cdk_descs, blued_descs
    # Split according to scheme
    if split == 'random':
        train_conns, test_conns = train_test_split(data.connectivity.unique(), test_size=0.2, random_state=SEED)
        train_conns, val_conns = train_test_split(train_conns, test_size=0.125, random_state=SEED)  # 0.125 x 0.8 = 0.1
        train_ids = np.where(data['connectivity'].isin(train_conns))[0]
        valid_ids = np.where(data['connectivity'].isin(val_conns))[0]
        test_ids = np.where(data['connectivity'].isin(test_conns))[0]
    elif split == 'scaffold':
        train_ids, valid_ids, test_ids = ScaffoldSplitter().train_valid_test_split(data, data[smilescolumn],
                                                                                   frac_train=0.70,
                                                                                   frac_valid=0.1,
                                                                                   frac_test=0.2,
                                                                                   random_state=SEED,
                                                                                   logfile=logfile)
    else:
        raise NotImplementedError(f'Split {split} is not implemented.')
    # Extract dependent variable
    X_train = descs.loc[descs.index.isin(train_ids)]
    X_valid = descs.loc[descs.index.isin(valid_ids)]
    X_test = descs.loc[descs.index.isin(test_ids)]
    y_train = data.loc[descs.index.isin(train_ids)][endpoint]
    y_valid = data.loc[descs.index.isin(valid_ids)][endpoint]
    y_test = data.loc[descs.index.isin(test_ids)][endpoint]
    # Scale the data
    if endpoint == 'logs':
        scaler = StandardScaler()
    elif endpoint == 'deltalogs':
        scaler = MaxAbsScaler()
    else:
        raise ValueError(f'Endpoint {endpoint} is not supported.')
    X_train.loc[:, :] = scaler.fit_transform(X_train)
    X_valid.loc[:, :] = scaler.transform(X_valid)
    X_test.loc[:, :] = scaler.transform(X_test)
    # Save scaler and dropped features for inference
    ml2json.to_json(scaler, f'{output}_scaler.json')
    with open(f'{output}_dropped_features.json', 'w') as oh:
        json.dump(X_train.loc[:, X_train.std() <= .1].columns.tolist(), oh)
    X_train = X_train.loc[:, X_train.std() > .1].astype('float64')
    X_valid = X_valid.loc[:, X_train.columns].astype('float64')
    X_test = X_test.loc[:, X_train.columns].astype('float64')
    del scaler
    # Train the model
    device = torch.device('cpu')
    if torch.cuda.is_available():
        device = torch.device('cuda')
    # Set fixed random number seed
    torch.manual_seed(SEED)
    last_save = 0  # last epoch at which the model was written to disk
    # Initialize the MLP
    mlp = MLPRegressor(X_train.shape[1]).to(device)
    print(mlp)
    print(f'Number of trainable parameters: {mlp.num_params:n}')
    # mlp = torch.compile(mlp)
    # Define the loss function and optimizer
    loss_function = nn.MSELoss()
    optimizer = torch.optim.Adam(mlp.parameters(), lr=LEARNING_RATE)
    best_loss = float('-inf')
    # Run the training loop
    for epoch in range(0, 100_000):
        train_inputs, train_targets = torch.from_numpy(X_train.values), torch.from_numpy(y_train.values)
        train_inputs, train_targets = train_inputs.float().to(device), train_targets.float()
        train_targets = train_targets.reshape((train_targets.shape[0], 1)).to(device)
        # Zero the gradients
        optimizer.zero_grad()
        # Perform forward pass
        train_outputs = mlp(train_inputs, dropout=True)
        # Compute loss
        loss = loss_function(train_outputs, train_targets)
        # Perform backward pass
        loss.backward()
        # Perform optimization
        optimizer.step()
        # Detach tensors
        train_targets = train_targets.cpu().detach().numpy().reshape((train_targets.shape[0],))
        train_outputs = train_outputs.cpu().detach().numpy().reshape((train_outputs.shape[0],))
        # Print training matrics
        print(f'Epoch {epoch + 1} (training):\tMSE Loss: {mean_squared_error(train_targets, train_outputs):.3f}\t'
              f'R²: {r2_score(train_targets, train_outputs):.3f}\tPearson\'s r: {pearsonr(train_targets, train_outputs).statistic:.3f}',
              file=logfile, flush=True)
        print(f'Epoch {epoch + 1} (training):\tMSE Loss: {mean_squared_error(train_targets, train_outputs):.3f}\t'
              f'R²: {r2_score(train_targets, train_outputs):.3f}\tPearson\'s r: {pearsonr(train_targets, train_outputs).statistic:.3f}')
        # Obtain and print validation metrics
        val_inputs, val_targets = torch.from_numpy(X_valid.values), torch.from_numpy(y_valid.values)
        val_inputs, val_targets = val_inputs.float(), val_targets.float()
        val_targets = val_targets.reshape((val_targets.shape[0], 1))
        val_outputs = mlp(val_inputs.to(device))
        val_targets = val_targets.cpu().detach().numpy().reshape((val_targets.shape[0],))
        val_outputs = val_outputs.cpu().detach().numpy().reshape((val_outputs.shape[0],))
        print(pd.Series(val_outputs).isnull().any())
        # Determine if converged
        # conv_criteria = mean_squared_error(val_targets, val_outputs)
        conv_criteria = pearsonr(val_targets, val_outputs).statistic
        print(f'\tEpoch {epoch + 1} (validation):\tMSE Loss: {mean_squared_error(val_targets, val_outputs):.3f}\t'
              f'R²: {r2_score(val_targets, val_outputs):.3f}\tPearson\'s r: {pearsonr(val_targets, val_outputs).statistic:.3f}',
              file=logfile, flush=True)
        print(f'\tEpoch {epoch + 1} (validation):\tMSE Loss: {mean_squared_error(val_targets, val_outputs):.3f}\t'
              f'R²: {r2_score(val_targets, val_outputs):.3f}\tPearson\'s r: {pearsonr(val_targets, val_outputs).statistic:.3f}')
        if conv_criteria > best_loss:
            torch.save(mlp.state_dict(), f'{output}_model.pkg')
            best_loss = conv_criteria
            last_save = epoch
        # Early stopping
        if (epoch >= MIN_EPOCHS) and (epoch - last_save > PATIENCE):
            print(
                f'No loss improvement in the last {PATIENCE} epochs: max Pearson\'s r = {best_loss:.6f} from epoch {last_save}',
                file=logfile, flush=True)
            print(
                f'No loss improvement in the last {PATIENCE} epochs: max Pearson\'s r = {best_loss:.6f} from epoch {last_save}')
            break
    stop_time = time.monotonic()
    # Process is complete.
    print(f'Training process has finished after {timedelta(seconds=stop_time - start_time)}',
          file=logfile, flush=True)
    print(f'Training process has finished after {timedelta(seconds=stop_time - start_time)}')
    # Reload model
    mlp.load_state_dict(torch.load(f'{output}_model.pkg'))
    test_inputs, test_targets = torch.from_numpy(X_test.values), torch.from_numpy(y_test.values)
    test_inputs, test_targets = test_inputs.float().to(device), test_targets.float()
    test_targets = test_targets.reshape((test_targets.shape[0], 1)).to(device)
    test_outputs = mlp(test_inputs.to(device), dropout=False)
    test_targets = test_targets.cpu().detach().numpy().reshape((test_targets.shape[0],))
    test_outputs = test_outputs.cpu().detach().numpy().reshape((test_outputs.shape[0],))
    # Print testing metrics
    print(f'Epoch {last_save} (test):\tMSE Loss: {mean_squared_error(test_targets, test_outputs):.3f}\t'
          f'R²: {r2_score(test_targets, test_outputs):.3f}\tPearson\'s r: {pearsonr(test_targets, test_outputs).statistic:.3f}',
          file=logfile, flush=True)
    print(f'Epoch {last_save} (test):\tMSE Loss: {mean_squared_error(test_targets, test_outputs):.3f}\t'
          f'R²: {r2_score(test_targets, test_outputs):.3f}\tPearson\'s r: {pearsonr(test_targets, test_outputs).statistic:.3f}')
    logfile.close()


@main.command(help='Convert a trained SolpHedge model from PyTorch to NumPy for inference.',
              context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input', type=str, required=True, multiple=False,
              default=None, nargs=1, show_default=False,
              help='File containing weights of the trained PyTorch model with extension.')
@click.option('-o', '--output', type=str, required=True, multiple=False,
              default=None, nargs=1, show_default=False,
              help='Name of the output pickle file (containing weights as NumPy arrays) without extension.')
def convert(input, output):
    """CLI to convert SolpHedge trained models for inference."""
    # Ensure input file exists
    if not os.path.exists(input):
        raise FileNotFoundError(f'File does not exist: {input}')
    out_path = output.rstrip('.pkg') + '.pkl'
    # Ensure output file does not exist or ask if overwriting
    if os.path.exists(out_path):
        while (confirm := raw_input('Output file already exists. Overwrite? (y/N)\n')).lower() not in ['y', 'n']:
            pass
        if confirm.lower() == 'n':
            return
    # Read weights of the PyTorch model
    weights = torch.load(input)
    out_weights = {key: mat.numpy(force=True) for key, mat in weights.items()}
    with open(out_path, 'wb') as oh:
        pickle.dump(out_weights, oh)
    print(f'Weights successfully exported to: {out_path}')


if __name__ == '__main__':
    main()