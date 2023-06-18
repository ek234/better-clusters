import os
import numpy as np
import pandas as pd
import pickle as pkl
from enum import Enum, auto
from sklearn.ensemble import AdaBoostClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, classification_report

root_dir = "./extracted/"
input_file = os.path.join(root_dir, "all_inputs.csv")
model_file = os.path.join(root_dir, "model.pkl")


class Event_types(str, Enum):
    hbbbar = "hbbbar"
    hlvjj = "hlvjj"
    ttbar_semilep = "ttbar_semilep"
    hgg = "hgg"
    ttbar_had = "ttbar_had"
    hjj = "hjj"
    ttbar_lep = "ttbar_lep"
    wjj = "wjj"
    zjj = "zjj"


content_df = pd.read_csv(input_file)

# pt cut
minpt = 10
content_df = content_df.loc[content_df["pt"] >= minpt]

labels = content_df[content_df.columns[0]].to_numpy()
features = content_df[content_df.columns[1:]].to_numpy().astype(np.float32)


labels_to_consider = [
    Event_types.hbbbar,
    Event_types.hlvjj,
    Event_types.ttbar_semilep,
    Event_types.hgg,
    Event_types.ttbar_had,
    Event_types.hjj,
    Event_types.ttbar_lep,
    Event_types.wjj,
    Event_types.zjj,
]
idx_to_consider = np.where(np.isin(labels, [lbl.value for lbl in labels_to_consider]))
features, labels = features[idx_to_consider], labels[idx_to_consider]

print(np.asanyarray(np.unique(labels, return_counts=True)).T)
features = np.nan_to_num(features, copy=False)

train_x, val_x, train_y, val_y = train_test_split(features, labels, test_size=0.3)

clf = AdaBoostClassifier(n_estimators=500)
clf.fit(train_x, train_y)
with open(model_file, "wb") as file:
    pkl.dump(clf, file)

print()

train_score = clf.score(train_x, train_y)
print(f"train score: {train_score}")
train_pred = clf.predict(train_x)
print("train conf matrix")
print(confusion_matrix(train_y, train_pred))
print("train classification report")
print(
    classification_report(
        train_y, train_pred, target_names=[lbl.name for lbl in labels_to_consider]
    )
)

val_score = clf.score(val_x, val_y)
print(f"val score: {val_score}")
val_pred = clf.predict(val_x)
print("val conf matrix")
print(confusion_matrix(val_y, val_pred))
print("val classification report")
print(
    classification_report(
        val_y, val_pred, target_names=[lbl.name for lbl in labels_to_consider]
    )
)
