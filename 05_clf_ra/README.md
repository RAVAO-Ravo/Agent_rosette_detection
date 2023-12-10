# Classifieur "clf_RA"

Ce dossier est dédié à la création d'un classifieur utilisant des séquences fictives générées à partir d'espèces modèles. Le processus global comprend la génération de données simulées, la création de jeux de données de test et de validation, la recherche du meilleur modèle, la recherche des meilleurs paramètres, et enfin, la sauvegarde du meilleur modèle avec les meilleurs paramètres.

## Structure du Projet

1. **`datasets/`**
   - Dossier contenant les données (sous forme vectorisée) utilisées pour la construction du modèle.

2. **`fastqs/`**
   - Dossier contenant les barcodes du 29 février 2023 utilisés pour la création des datasets de test et de validation.

3. **`species/`**
   - Dossier contenant les fichiers FASTA des espèces modèles nécessaires à la génération de données simulées pour l'entraînement du modèle.

4. **`utilitaires.py`**
   - Fichier contenant des variables globales et des fonctions utilitaires nécessaires pour le projet.

5. **`make_realDataset.py`**
   - Fichier permettant la création du jeu de données de test et de validation à partir des barcodes 06 et 12 du dossier `fastqs`, suivi de la vectorisation dans le dossier `datasets`.

6. **`make_fakeDataset.py`**
   - Fichier permettant la création de séquences fictives avec des erreurs à partir des espèces du dossier `species`, suivi de la vectorisation dans le dossier `datasets`.

7. **`findBestModel.ipynb`**
   - Notebook destiné à la recherche du meilleur modèle à utiliser pour le classifieur.

8. **`findBestParams.ipynb`**
   - Notebook destiné à trouver les meilleurs paramètres pour le meilleur modèle identifié précédemment.

9. **`makeCLF_RA.py`**
   - Fichier permettant de sauvegarder le meilleur modèle avec les meilleurs paramètres sous forme d'un fichier pickle.

## Déroulement du Projet

1. Génération de séquences fictives :
   - Utilisation du fichier `make_fakeDataset.py` pour créer des séquences fictives à partir des espèces modèles du dossier `species` avec introduction d'erreurs. Les données sont ensuite vectorisées dans le dossier `datasets`.

2. Création des données de test et de validation :
   - Utilisation du fichier `make_realDataset.py` pour créer les données de test et de validation à partir des barcodes 06 et 12 du dossier `fastqs`. Les données sont vectorisées dans le dossier `datasets`.

3. Recherche du meilleur modèle :
   - Exploration des modèles à l'aide du notebook `findBestModel.ipynb` pour identifier le modèle le plus performant pour la tâche de classification.

4. Recherche des meilleurs paramètres :
   - Utilisation du notebook `findBestParams.ipynb` pour trouver les paramètres optimaux pour le meilleur modèle identifié précédemment.

5. Sauvegarde du meilleur modèle :
   - Utilisation du fichier `makeCLF_RA.py` pour sauvegarder le meilleur modèle avec les meilleurs paramètres dans un fichier pickle.

## Note Importante

Pour que le projet fonctionne correctement, assurez-vous d'inclure les fichiers `.so` issus de la compilation des `.pyx` du dossier "00_reads".