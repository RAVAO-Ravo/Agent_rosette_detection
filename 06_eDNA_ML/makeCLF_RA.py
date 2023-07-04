#!/bin/python3
#-*- coding:utf-8 -*-


###### Importations ######

import pandas as pd
import pickle as pk
from xgboost import XGBClassifier

###### Importations ######


###### Variables globales ######

encoder: dict = {"Anurofeca_richardsi": 0,
                 "Dermocystidium_salmonis": 1,
                 "Ichthyophonus_hoferi": 2,
                 "Pseudoperkinsus_tapetis": 3,
                 "Psorospermium_haeckelii": 4,
                 "Rhinosporidium_cygnus": 5,
                 "Rhinosporidium_seeberi": 6,
                 "Sphaeroforma_spB7": 7,
                 "Sphaeroforma_spCRG3": 8,
                 "Sphaerothecum_destruens": 9,
                 "unidentified": 10}
random_state: int = 42
species_amplicons: pd.DataFrame = pd.read_csv(filepath_or_buffer="./dataset/fakes_species_amplicons.csv")
randoms_amplicons: pd.DataFrame = pd.read_csv(filepath_or_buffer="./dataset/fakes_unidentified_amplicons.csv")
best_size: int = 44000
best_params: list = [400, 3, 0.2, 0.6, 0.6, 0.9]
model_name: str = "clf_RA.sav"

###### Variables globales ######


if __name__ == "__main__" :

	print("Chargement des données")
	data = pd.concat(objs=[species_amplicons, randoms_amplicons], axis=0).sample(frac=1.0, replace=False, random_state=random_state).reset_index(drop=True)
	x_train, y_train = data.iloc[:, :-1], data.iloc[:, -1].apply(lambda x: encoder[x])

	print("Entrainement du modèle")
	xgb = XGBClassifier(n_estimators=best_params[0],
						max_depth=best_params[1],
						learning_rate=best_params[2],
						subsample=best_params[3],
						colsample_bytree=best_params[4],
						colsample_bylevel=best_params[5], 
						random_state=random_state).fit(X=x_train.iloc[:best_size], y=y_train.iloc[:best_size])

	print("Sauvegarde du modèle")
	pk.dump(obj=xgb, file=open(file=model_name, mode="wb"))