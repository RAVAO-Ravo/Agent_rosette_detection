# Reads

Il s'agit d'un ensemble de modules Cython conçus pour la manipulation efficace de fichiers au format FASTQ. Le module principal, `reads.pyx`, dépend de plusieurs éléments pour assurer son bon fonctionnement.

## Contenu du Dossier

1. **`reads.pyx`**
   - Ce module central permet la manipulation des fichiers FASTQ. Il utilise les fonctionnalités fournies par les autres fichiers Cython.

2. **`config.pyx`**
   - Fichier Cython contenant des éléments de configuration essentiels au bon fonctionnement du module `reads.py`.

3. **`ahocorasick.pyx`**
   - Implémentation de l'algorithme Aho-Corasick, un algorithme efficace pour la recherche de motifs dans un texte. Ce fichier est essentiel pour les opérations de recherche dans les données FASTQ.

4. **`read.pyx`**
   - Implémentation d'une structure de données spécifique qui est utilisée par le module `reads.pyx`.

## Utilisation

1. Compiler les fichiers avec la commande `make` dans un terminal :
   
   ```shell
   make
   ```

2. Importez le module `reads` dans votre code Python :
   
   ```python
   from reads import Reads
   ```