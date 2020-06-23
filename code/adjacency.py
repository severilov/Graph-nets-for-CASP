import os
import numpy as np
import datetime
from tqdm import tqdm
from biopandas.pdb import PandasPdb
from Bio.PDB import NeighborSearch, PDBParser

ppdb = PandasPdb()
SEARCH_RADIUS = 6
KMIN_DISTANCE_BETWEEN_ATOMS = 0.01
WATER_RESIDUES = {'HOH', 'HHO', 'OHH', 'H2O', 'OH2', 'WAT', 'TIP',
                  'TIP3', 'TIP4', 'TIP3P', 'TIP4P', 'SOL'}

def get_elements_radii(elements_radii_path):
    return dict(map(
        lambda l: (l.strip().split(': ')[0], float(l.strip().split(': ')[1])),
        open(elements_radii_path).readlines()))


def get_chain_id(atom):
    return atom.get_parent().get_parent()._id


def get_residue_number(atom):
    return atom.get_parent()._id[1]


def get_residue_name(atom):
    return atom.get_parent().resname


def has_water_residue(atom):
    return get_residue_name(atom) in WATER_RESIDUES


def is_hydrogen(atom):
    return atom.element == 'H'


def have_covalent_bond(first, second, elements_radii):
    distance_2 = np.linalg.norm(first.get_coord() - second.get_coord(), ord=2) ** 2

    if has_water_residue(first) and has_water_residue(second):
        return distance_2 < 1.21

    # then, we don't connect water with anything else
    if has_water_residue(first) ^ has_water_residue(second):
        return False

    # another rule, we don't connect things from different chains
    if get_chain_id(first) != get_chain_id(second):
        return False

    if is_hydrogen(first) and is_hydrogen(second):
        return False

    if (is_hydrogen(first) or is_hydrogen(second)) and distance_2 >= 1.21:
        return False

    # we don't create S-S bridges
    if first.name.startswith('SG') and second.name.startswith('SG'):
        return False

    # we don't connect things that are far in the sequence
    if abs(get_residue_number(first) - get_residue_number(second)) > 1:
        return False

    # here we don't want to have a connection between two "H" in a water
    # but there are some forcefields where it can be possible

    distance_max_2 = (elements_radii.get(first.element, -1) + elements_radii.get(second.element, -1)) * 0.6
    distance_max_2 *= distance_max_2

    # we don't consider pairs which are too close also because that could be a mistake
    return (distance_2 <= distance_max_2) and (distance_2 > KMIN_DISTANCE_BETWEEN_ATOMS)

def AdjacencyMatrix(models_path, target_name, model_name):

    pdb_file_path = models_path + target_name + '/' + model_name
    temp = ppdb.read_pdb(pdb_file_path)
    protein = ppdb.df['ATOM']
    # some model files have empty data about molecule (e.g. MIG_FROST_AL1 in CASP7)
    if len(protein) == 0:
        print('BAD ' + pdb_file_path)
        return pdb_file_path
    # it may get error because some model files have nans instead of coordinates (e.g. in CASP8)
    # some model files have HEATOM, it may cause error too
    else:
        #tempor = protein[['atom_number','atom_name', 'residue_name', 'residue_number']]
        tempor = protein[['atom_name', 'residue_name', 'residue_number']].T.to_dict('list')
        keys=list(tempor.keys())
        values=list(tempor.values())


        #atom_to_pos = tempor.set_index(['atom_name', 'residue_name', 'residue_number']).to_dict('index')

        if not os.path.exists(models_path + target_name + '/adjacency'):
            os.mkdir(models_path + target_name + '/adjacency')

        f = open(models_path + target_name + '/adjacency/' + model_name + '_adj.txt', "w")

        elements_radii_path = './elements_radii.txt'
        elements_radii = get_elements_radii(elements_radii_path)
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        structure = parser.get_structure(model_name, pdb_file_path)
        atoms = [atom for atom in structure.get_atoms()]
        neighbors_searcher = NeighborSearch(atoms)
        pairs = neighbors_searcher.search_all(radius=SEARCH_RADIUS)
        #covalent_bonds = []
        for atom_1, atom_2 in pairs:
            if have_covalent_bond(atom_1, atom_2, elements_radii):
                atom_id1 = keys[values.index([atom_1.name,
                                          get_residue_name(atom_1),
                                          get_residue_number(atom_1)])]
                atom_id2 = keys[values.index([atom_2.name,
                                          get_residue_name(atom_2),
                                          get_residue_number(atom_2)])]
                '''
                # если хотим 'atom_number', а не позицию
                atom_id1 = atom_to_pos[atom_1.name,
                                        get_residue_name(atom_1),
                                        get_residue_number(atom_1)]['atom_number']
                atom_id2 = atom_to_pos[atom_2.name,
                                        get_residue_name(atom_2),
                                        get_residue_number(atom_2)]['atom_number']
                '''
                f.write(str(atom_id1) + '_' + str(atom_id2) + ' ')
                # возможно стоит записывать еще прочие признаки, как residue_number, distance, atom_name, residue_name
            #f.write('\n')
        f.close()
        return 'OK'

def MakeAdjForAll(models_path, targets_path, casp_name):
    print('/'*100 + '\nStarting ' + casp_name + '\n' + '/'*100)
    models_count = 0
    unwanted = {'.DS_Store', 'adjacency_old', 'adjacency'}
    targets = os.listdir(models_path)
    target_names = sorted([s.replace('.pdb', '') for s in targets])[1:] # [1:] т.к. считывает системный файл .DS_store
    targets = sorted(target_names)

    bad_models = []
    for i, target_name in enumerate(targets[114:]):
        models = sorted(os.listdir(models_path + target_name))
        models = [model for model in models if model not in unwanted]
        models_count += len(models)
        print('Target [{}/{}] {} started at {}'.format(114+i+1, len(targets),
                                                       target_name,
                                                       datetime.datetime.now().time()))
        for model in tqdm(models):
            bad = AdjacencyMatrix(models_path, target_name, model)
            if bad != 'OK':
                bad_models.append(bad)

    print('/'*100)
    print('TOTAL for {}\ntargets: {}, models: {}'.format(casp_name,
                                                        len(targets),
                                                        models_count))
    return bad_models


if __name__ == "__main__":
    # for colab
    #models_path = '/content/drive/My Drive/CASP/CASP11/models/'
    #target_path = '/content/drive/My Drive/CASP/CASP11/targets/'

    models_paths = ['../../data/CASP7/models/', '../../data/CASP8/models/',
                    '../../data/CASP9/models/', '../../data/CASP10/models/',
                    '../../data/CASP11/models/','../../data/CASP12/models/']
    targets_paths = ['../../data/CASP7/targets/','../../data/CASP8/targets/',
                     '../../data/CASP9/targets/', '../../data/CASP10/targets/',
                     '../../data/CASP11/targets/','../../data/CASP12/targets/']
    competitions = ['CASP7', 'CASP8', 'CASP9', 'CASP10', 'CASP11', 'CASP12']

    # done for CASP7, CASP8, CASP10, CASP11, CASP12
    for i in range(len(competitions)):
        bad_models = MakeAdjForAll(models_paths[i], targets_paths[i], competitions[i])
        print('!!! BAD MODELS:{}'.format(bad_models))
        print('/'*100)
