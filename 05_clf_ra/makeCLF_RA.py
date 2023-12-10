#!/bin/python3
#-*- coding:utf-8 -*-

###### Importations ######

# Importation des modules
import pandas as pd
import pickle as pk
from typing import List, Union

# Ignorer les avertissements spécifiques de la bibliothèque XGBoost
import warnings
warnings.filterwarnings(action="ignore", module="xgboost")
import xgboost as xgb_

# Importation des éléments de "utilitaire.py"
from utilitaire import (
	SPECIES_AMPLICONS_FILE,
	RANDOMS_AMPLICONS_FILE,
	CLF_FILE,
	RANDOM_STATE,
	encode_sp
)

###### Importations ######

###### Variables globales ######

SPECIES_AMPLICONS: pd.DataFrame = pd.read_csv(filepath_or_buffer=SPECIES_AMPLICONS_FILE)
RANDOMS_AMPLICONS: pd.DataFrame = pd.read_csv(filepath_or_buffer=RANDOMS_AMPLICONS_FILE)
BEST_SIZE: int = 88000
BEST_PARAMS: List[Union[int, float]] = [800, 20, 0.3, 0.30, 0.70, 1.0]

###### Variables globales ######

if __name__ == "__main__" :
	print("Chargement des données")
	data = pd.concat(objs=[SPECIES_AMPLICONS, RANDOMS_AMPLICONS], axis=0).sample(frac=1.0, replace=False, random_state=RANDOM_STATE).reset_index(drop=True)
	x_train, y_train = data.iloc[:, :-1], data.iloc[:, -1].apply(encode_sp)

	print("Entrainement du modèle")
	xgb = xgb_.XGBClassifier(n_estimators=BEST_PARAMS[0],
							 max_depth=BEST_PARAMS[1],
							 learning_rate=BEST_PARAMS[2],
							 subsample=BEST_PARAMS[3],
							 colsample_bytree=BEST_PARAMS[4],
							 colsample_bylevel=BEST_PARAMS[5], 
							 random_state=RANDOM_STATE).fit(X=x_train.iloc[:BEST_SIZE], y=y_train.iloc[:BEST_SIZE])

	print("Sauvegarde du modèle")
	pk.dump(obj=xgb, file=open(file=CLF_FILE, mode="wb"))