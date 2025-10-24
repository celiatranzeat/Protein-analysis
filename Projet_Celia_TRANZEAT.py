"""
Module contenant les fonctions et les tests des fonctions 
"""

##########################################################################
###IMPORTATION DES MODULES:
##########################################################################

import Bio
from Bio import ExPASy
from Bio import SwissProt
from Bio.SeqUtils.ProtParam import ProteinAnalysis, ProtParamData
import matplotlib.pyplot as plt
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
import unittest


##########################################################################
##FONCTIONS:
##########################################################################

def fetch_uniprot_record(protein_id):
    """
    Récupère un enregistrement UniProt à partir d'un identifiant ExPASy.

    Parameters
    ----------
    protein_id : str
        Identifiant de la protéine dans la base ExPASy

    Returns
    -------
    record: SeqRecord
        Objet SeqRecord contenant la séquence et les annotations associées à la protéine.
    """
    handle = ExPASy.get_sprot_raw(protein_id)
    record = SwissProt.read(handle)
    handle.close()
    return record

def get_protein_dictinfo(record):
    """
    Extrait les informations principales de la protéine
    et les retourne sous forme de dictionnaire.

    Parameters
    ----------
    record : SeqRecord
        Objet SeqRecord contenant la séquence et les annotations de la protéine.

    Returns
    -------
    dict
        Dictionnaire avec les clés suivantes :

        - "Uniprot_id" : Identifiant UniProt de la protéine.
        - "Gene_name" : Nom du gène associé.
        - "Accessions" : Liste des numéros d'accession.
        - "Sequence Length" : Longueur de la séquence protéique.
        - "Description" : Description de la protéine.
        - "Sequence" : Séquence protéique en tant que chaîne de caractères.
        - "Taxonomy" : Taxonomie associée.
        - "Comments" : Informations sur la localisation subcellulaire.
    """
    dict_info = {
        "Uniprot_id": record.entry_name,
        "Gene_name": record.gene_name,
        "Accessions": record.accessions,
        "Sequence Length": record.sequence_length,
        "Description": record.description,
        "Sequence": record.sequence,
        "Taxonomy": record.organism_classification,
        "Comments subcellular location" : record.comments[0]}

    return dict_info



def display_dictinfo(dict_info):
    """
    Affiche le contenu d'un dictionnaire d'informations sous forme lisible.

    Format d’affichage
    ------------------
    Pour chaque paire (clé, valeur), la fonction imprime :
        clé
        valeur
    suivies d’un retour à la ligne.

    Parameters
    ----------
    dict_info : dict
        Dictionnaire d’informations

    Returns
    -------
    None
        La fonction imprime simplement le texte.
    """
    info = ""
    for key, value in dict_info.items():
        info += f"{key}\n{value}\n"
    print(info)


def get_aa_composition(record):
    """
    Calcule la composition en acides aminés d'une protéine (fraction par acide aminé).

    Parameters
    ----------
    record : SeqRecord
        Objet SeqRecord contenant la séquence et les annotations de la protéine.

    Returns
    -------
    dict[str: float]
        Dictionnaire {AA_1lettre: fraction} où chaque valeur est comprise entre 0.0 et 1.0.
    """
    X = ProteinAnalysis(record.sequence)
    composition = X.get_amino_acids_percent()
    return composition


def display_composition_bar(composition):
    """
    Affiche sous forme d'histogramme le pourcentage des acides aminés.

    Parameters
    ----------
    composition : dict[str: float]
        Dictionnaire {AA_1lettre: fraction} où chaque valeur est comprise entre 0.0 et 1.0.
        Les valeurs sont multipliées par 100 pour l’affichage.

    Returns
    -----
    Affiche l'histogramme matplotlib.
    """
    names = list(composition.keys())
    values = [v * 100 for v in composition.values()]
    plt.figure(figsize=(8, 5))
    plt.bar(names, values)
    plt.title("Composition en acides aminés")
    plt.xlabel("Acides aminés")
    plt.ylabel("Pourcentage")
    plt.show()

def get_hydrophobicity(record):
    """
    Calcule le profil d’hydrophobicité d’une protéine.

    Parameters
    ----------
    record : SeqRecord
        Objet SeqRecord contenant la séquence et les annotations de la protéine.

    Returns
    -------
    list[float]
        Profil d’hydrophobicité (longueur ≈ len(sequence) - window + 1).

    Interprétation
    --------------
    - Valeurs positives: segments plutôt hydrophobes.
    - Valeurs négatives : segments plutôt hydrophiles.
    """
    X = ProteinAnalysis(record.sequence)
    hydrophobicity = X.protein_scale(param_dict=ProtParamData.kd, window=9)
    return hydrophobicity


def display_hydrophobicity(hydrophobicity):
    """
    Trace un profil d’hydrophobicité de la protéine en fonction de
    la position des acides aminés.

    Parameters
    ----------
    hydrophobicity : list[float]
        Liste de valeurs d’hydrophobicité à tracer.

    Returns
    -----
    Affiche une figure "plot" matplotlib avec :
      - abscisse : positions des acides aminés,
      - ordonnée : valeurs d’hydrophobicité.
    """
    position = list(range(1, len(hydrophobicity) + 1))
    plt.plot(position, hydrophobicity)
    plt.xlabel("Position de l'acide aminé")
    plt.ylabel("Score d'hydrophobicité")
    plt.title("Profil d'hydrophobicité")
    plt.show()

def ph_charge(record):
    """
    Calcule la courbe « charge nette en fonction du pH » d'une protéine.

    Le calcul s'appuie sur Bio.SeqUtils.IsoelectricPoint.IsoelectricPoint et
    évalue la charge nette aux pH 0.0 jusqu'à 14.0 par pas de 0.5.

    Le point isoélectrique (pI) est approximativement le pH où la charge croise 0

    Parameters
    ----------
    record : SeqRecord
        Objet SeqRecord contenant la séquence et les annotations de la protéine.

    Returns
    -------
    tuple (pH, charge)
        pH : list[float]
        charge : list[float]
    """
    seq = record.sequence
    iep = IsoelectricPoint(seq)
    pH = [x / 2 for x in range(0, 29)]
    charge = [iep.charge_at_pH(x) for x in pH]
    return (pH, charge)


def display_ph_charge(list_values):
    pH, charge = list_values
    plt.plot(pH, charge)
    plt.xlabel("pH")
    plt.ylabel("Charge nette")
    plt.show()

def find_nglyco_motif(record):
    """
    Repère les motifs potentiels de N-glycosylation dans une séquence protéique.

    Définition du motif
    -------------------
    Un motif de N-glycosilation est: N-X-[S/T] avec X ≠ P

    Parameters
    ----------
    record : SeqRecord
        Objet SeqRecord contenant la séquence et les annotations de la protéine.

    Returns
    -------
    list[int]
        Liste des positions du résidu 'N' pour chaque motif détecté.
    """
    seq = str(record.sequence)
    pos = []

    for i in range(len(seq) - 2):
        a, b, c = seq[i], seq[i + 1], seq[i + 2]

        if a == "N" and b != "P" and c in ("S", "T"):
            pos.append(i + 1)

    return pos


##########################################################################
##TESTS DES FONCTIONS
##########################################################################

protein_id = "O23729" 
record = fetch_uniprot_record(protein_id) 

class TestFetchUniprotRecord(unittest.TestCase):
    """
    Tests unitaires pour la fonction fetch_uniprot_record().
    Vérifie que l'appel avec un identifiant UniProt valide retourne bien
    un objet SwissProt.Record correspondant à la bonne protéine.
    """
    def test_valid_id(self):
        self.assertIsInstance(record, SwissProt.Record)
        self.assertEqual(record.accessions[0], protein_id)


class TestGetProteinDictInfo(unittest.TestCase):
    """
    Tests unitaires pour la fonction get_protein_dictinfo().
    Vérifie que le dictionnaire renvoyé contient toutes les informations attendues,
    avec les bons types et une cohérence avec l'objet record source.
    """
    def setUp(self):
        """
        Initialise les données de test à partir du record.
        """
        self.info = get_protein_dictinfo(record)

    def test_keys_presence(self):
        """
        Vérifie que toutes les clés attendues sont présentes dans le dictionnaire.
        """
        expected_keys = [
            "Uniprot_id", "Sequence", "Sequence Length", "Gene_name",
            "Accessions", "Description", "Taxonomy"
        ]
        for key in expected_keys:
            self.assertIn(key, self.info)

    def test_sequence_length(self):
        """
        Vérifie que la longueur indiquée correspond bien à la longueur réelle de la séquence.
        """
        self.assertEqual(self.info["Sequence Length"], len(self.info["Sequence"]))

    def test_field_types(self):
        """
        Vérifie le type de chaque valeur du dictionnaire
        """
        self.assertIsInstance(self.info, dict)
        self.assertIsInstance(self.info["Uniprot_id"], str)
        self.assertIsInstance(self.info["Accessions"], list)
        self.assertIsInstance(self.info["Sequence Length"], int)
        self.assertIsInstance(self.info["Description"], str)
        self.assertIsInstance(self.info["Sequence"], str)
        self.assertIsInstance(self.info["Taxonomy"], list)

    def test_values_match_record(self):
        """
        Vérifie que certaines valeurs du dictionnaire correspondent au record source.
        """
        self.assertEqual(self.info["Accessions"], record.accessions)
        self.assertEqual(self.info["Uniprot_id"], record.entry_name)
        

class TestGetAAComposition(unittest.TestCase):
    """
    Tests unitaires pour la fonction get_aa_composition().
    Vérifie que la composition en acides aminés est correcte en type, format et cohérence.
    """

    def setUp(self):
        """
        Initialise la composition en acides aminés à partir du record.
        """
        self.composition = get_aa_composition(record)

    def test_return_type(self):
        """
        Vérifie que le résultat est un dictionnaire.
        """
        self.assertIsInstance(self.composition, dict)

    def test_keys_are_amino_acids(self):
        """
        Vérifie que chaque clé est une lettre d’acide aminé.
        """
        for aa in self.composition.keys():
            self.assertIsInstance(aa, str)
            self.assertEqual(len(aa), 1)

    def test_values_are_floats_between_0_and_1(self):
        """
        Vérifie que chaque valeur est un flottant compris entre 0 et 1.
        """
        for val in self.composition.values():
            self.assertIsInstance(val, float)
            self.assertGreaterEqual(val, 0.0)
            self.assertLessEqual(val, 1.0)

class TestGetHydrophobicity(unittest.TestCase):
    """
    Tests unitaires pour la fonction get_hydrophobicity().
    Vérifie le type, la longueur et la validité des valeurs calculées de l’hydrophobicité.
    """

    def setUp(self):
        """
        Initialise les valeurs d’hydrophobicité pour le record.
        """
        self.window = 9
        self.h = get_hydrophobicity(record)

    def test_return_type(self):
        """
        Vérifie que le résultat est une liste.
        """
        self.assertIsInstance(self.h, list)

    def test_length(self):
        """
        Vérifie que la longueur correspond à la séquence moins la fenêtre + 1.
        """
        expected_length = len(record.sequence) - self.window + 1
        self.assertEqual(len(self.h), expected_length)
            

class TestPHCharge(unittest.TestCase):
    """
    Tests unitaires pour la fonction ph_charge().
    Vérifie que les listes de pH et de charge retournées sont cohérentes en type,
    en longueur et en contenu numérique.
    """
    def test_ph_charge(self):
        pH, charge = ph_charge(record)

        self.assertIsInstance(pH, list)
        self.assertIsInstance(charge, list)
        self.assertEqual(len(pH), 29)              
        self.assertEqual(len(charge), len(pH))

        for v in charge:
            self.assertIsInstance(v, float)


class TestFindNGlycoMotif(unittest.TestCase):
    """
    Tests unitaires pour la fonction find_nglyco_motif().
    Vérifie que la liste des positions retournées est correcte en type
    et que chaque triplet détecté correspond bien au motif N-X-[S/T] attendu.
    """

    def test_return_type(self):
        """
        Vérifie que le résultat est une liste d’entiers.
        """
        seq = str(record.sequence)
        positions = find_nglyco_motif(record)
        self.assertIsInstance(positions, list)
        for p in positions:
            self.assertIsInstance(p, int)

    def test_motif_structure(self):
        """
        Vérifie que chaque triplet détecté correspond au motif N-X-[S/T].
        """
        seq = str(record.sequence)
        positions = find_nglyco_motif(record)
        for p in positions:
            tri = seq[p-1 : p+2]
            self.assertEqual(len(tri), 3)
            self.assertEqual(tri[0], "N")
            self.assertNotEqual(tri[1], "P")
            self.assertIn(tri[2], ("S", "T"))

            
def run_tests():
    """
    Fonction qui exécute tous les tests
    """
    test_suite = unittest.TestSuite()
    loader = unittest.TestLoader()

    test_suite.addTests(loader.loadTestsFromTestCase(TestFetchUniprotRecord))
    test_suite.addTests(loader.loadTestsFromTestCase(TestGetProteinDictInfo))
    test_suite.addTests(loader.loadTestsFromTestCase(TestGetAAComposition))
    test_suite.addTests(loader.loadTestsFromTestCase(TestGetHydrophobicity))
    test_suite.addTests(loader.loadTestsFromTestCase(TestPHCharge))
    test_suite.addTests(loader.loadTestsFromTestCase(TestFindNGlycoMotif))

    test_runner = unittest.TextTestRunner(verbosity=2)
    test_runner.run(test_suite)



