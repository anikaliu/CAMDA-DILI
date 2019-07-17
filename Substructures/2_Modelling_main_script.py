import numpy as np
import pandas as pd
import sklearn as sk
from sklearn import svm
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
import sys
import csv
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import SelectFromModel


def import_dataset():
    with open('binvec_ppv0.csv', newline='') as csvfile_in:
        bin_vec = pd.read_csv(csvfile_in, sep=',')
    print("Dataset: ", bin_vec.head())
    return bin_vec
    print("Dataset Length: ", len(bin_vec))
    print("Dataset Shape: ", bin_vec.shape)


def extract_x_y(bin_vec):
    x = bin_vec.values[:, 2:]
    y = bin_vec.values[:, 0]
    y = y.astype('int')
    return (x, y)


def split_dataset(x, y):
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3)
    print(x, y, X_train, X_test, y_train, y_test)
    return x_train, x_test, y_train, y_test


def model_assessment(model_name, model, x, y):
    x_train, x_test, y_train, y_test = split_dataset(x, y)
    model.fit(x_train, y_train)

    predicted = model.predict(x_test)
    conf_matrix = confusion_matrix(y_test, predicted)
    print(model_name, "conf_matrix=\n", conf_matrix)

    cv_score = cross_val_score(model, x, y, cv=5)
    print(model_name, "cv_score = ", cv_score)


def model_grid_search(model_name, model, x, y, param_grid):
    grid_search = GridSearchCV(model, param_grid, cv=StratifiedKFold(n_splits=5))
    grid_search.fit(x, y)
    print(model_name, "with Grid search best params", grid_search.best_params_)

    x_train, x_test, y_train, y_test = split_dataset(x, y)
    predicted = grid_search.predict(x_test)
    conf_matrix = confusion_matrix(y_test, predicted)
    print(model_name, "conf_matrix=\n", conf_matrix)

    # Now we know which parameters are best for the model.
    # We can use those to create a new model of the same kind.
    # For example:
    # m = DecisionTreeClassifier()
    # best_params = model_grid_search(".....", m, x, y, param_grid)
    # m2 = DecisionTreeClassifier(**best_params)

    return grid_search.best_params_


def kfold(model, X, Y):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)

    acc = []
    bal_acc = []
    rec_p = []
    rec_n = []
    pre_p = []
    pre_n = []
    f1_p = []
    f1_n = []
    roc_auc = []

    for train, test in skf.split(X, Y):
        model.fit(X[train], Y[train])
        acc.append(accuracy_score(Y[test], model.predict(X[test])))
        bal_acc.append(balanced_accuracy_score(Y[test], model.predict(X[test])))
        rec_p.append(recall_score(Y[test], model.predict(X[test]), pos_label=1))
        rec_n.append(recall_score(Y[test], model.predict(X[test]), pos_label=0))
        pre_p.append(precision_score(Y[test], model.predict(X[test]), pos_label=1))
        pre_n.append(precision_score(Y[test], model.predict(X[test]), pos_label=0))
        f1_p.append(f1_score(Y[test], model.predict(X[test]), pos_label=1))
        f1_n.append(f1_score(Y[test], model.predict(X[test]), pos_label=0))
        roc_auc.append(roc_auc_score(Y[test], model.predict(X[test])))

    print('accuracy: ', np.mean(acc))
    print('balanced accuracy: ', np.mean(bal_acc))
    print('recall positives: ', np.mean(rec_p))
    print('recall negatives: ', np.mean(rec_n))
    print('precision positives: ', np.mean(pre_p))
    print('precision negatives: ', np.mean(pre_n))
    print('f1 positives: ', np.mean(f1_p))
    print('f1 negatives: ', np.mean(f1_n))
    print('ROC AUC: ', np.mean(roc_auc))

# Main programme----------------------------------------------------------------


bin_vec = import_dataset()
x, y = extract_x_y(bin_vec)

# 5 models with cross validation and no grid_search ---------------------------
print("--------------------DECISION TREE - GINI (CV)--------------------")
clf_gini = DecisionTreeClassifier(criterion="gini",
                                  random_state=100, max_depth=3, min_samples_leaf=5)
model_assessment("Gini", clf_gini, x, y)
kfold(clf_gini, x, y)

print("--------------------DECISION TREE - ENTROPY (CV)--------------------")
clf_entropy = DecisionTreeClassifier(criterion="entropy",
                                     random_state=100, max_depth=3, min_samples_leaf=5)
model_assessment("Entropy", clf_entropy, x, y)
kfold(clf_entropy, x, y)

print("-------------------- RANDOM FOREST (CV)-----------------------------")
clf_rf = RandomForestClassifier(max_depth=9, oob_score=True,
                                n_estimators=500, criterion='gini', random_state=2, min_samples_leaf=1,
                                class_weight='balanced_subsample', bootstrap=True, max_features='auto')
model_assessment("Random Forest", clf_rf, x, y)
kfold(clf_rf, x, y)
print("--------------------SVM (CV)----------------------------------------")
clf_svm = svm.SVC(gamma='scale')
model_assessment("SVM", clf_svm, x, y)
kfold(clf_svm, x, y)

print("-------------------XGBOOST (CV)------------------------------------")
clf_xgbc = XGBClassifier()
model_assessment("XGBoost", clf_xgbc, x, y)
kfold(clf_xgbc, x, y)

print("--------------------------------MODELS WITH CROSS VALIDATION and GRID SEARCH--------------------------------")

print("--------------------DECISION TREE - GINI (CV,GD)--------------------")
clf_gini2 = DecisionTreeClassifier(criterion="gini", random_state=100)
gini2_param_grid = {
    'max_depth': [3, 5, 7, 9],
    'min_samples_leaf': [3, 5, 7, 9],
    'max_features': ['auto', 'sqrt'],
}
gini_best_params = model_grid_search("Gini", clf_gini2, x, y, gini2_param_grid)
# Unpack the best params dict with **gini_best_params
clf_gini3 = DecisionTreeClassifier(criterion="gini",
                                   random_state=100, **gini_best_params)
kfold(clf_gini3, x, y)
model_assessment("Gini", clf_gini3, x, y)
print(gini_best_params)
print("--------------------DECISION TREE - ENTROPY (CV,GD)--------------------")
clf_entropy2 = DecisionTreeClassifier(criterion="entropy", random_state=100)
entropy2_param_grid = {
    'max_depth': [3, 5, 7, 9],
    'min_samples_leaf': [3, 5, 7, 9],
    'max_features': ['auto', 'sqrt'],
}
entropy_best_params = model_grid_search("entropy", clf_entropy2, x, y, entropy2_param_grid)
clf_entropy3 = DecisionTreeClassifier(criterion="entropy",
                                      random_state=100, **entropy_best_params)
kfold(clf_entropy3, x, y)
model_assessment("Entropy", clf_entropy3, x, y)
print(entropy_best_params)

print("--------------------RANDOM FOREST (CV,GD)------------------------------")
clf_rf2 = RandomForestClassifier(oob_score=True, criterion="gini",
                                 random_state=2, class_weight='balanced_subsample', bootstrap=True)
rf2_param_grid = {
    'max_depth': [3, 6, 9, 10, 12, 15, None],
    'max_features': ['auto', 'sqrt'],
    'min_samples_leaf': [1, 2, 3, 4],
    'min_samples_split': [2, 5, 10],
    'n_estimators': [20, 60, 100, 200, 300, 400, 500, 750, 1000]
}
rf2_best_params = model_grid_search("gini", clf_rf2, x, y, rf2_param_grid)
clf_rf3 = RandomForestClassifier(criterion="gini", random_state=100, **rf2_best_params)
kfold(clf_rf3, x, y)
model_assessment("Random Forest", clf_rf3, x, y)
print(rf2_best_params)

clf_rf3.fit(x, y)
importances = clf_rf3.feature_importances_
std = np.std([tree.feature_importances_ for tree in clf_rf3.estimators_],
             axis=0)
indices = np.argsort(importances)[::-1]

# Print the feature ranking
print("Feature ranking:")

for f in range(x.shape[1]):
    print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

# Plot the feature importances of the forest
plt.figure()
plt.title("Feature importances")
plt.bar(range(x.shape[1]), importances[indices],
        color="r", yerr=std[indices], align="center")
plt.xticks(range(x.shape[1]), indices)
plt.xlim([-1, x.shape[1]])
plt.show()

#df = pd.DataFrame(np.random.rand(10, 5), columns=['A', 'B', 'C', 'D', 'E'])
# df.plot.box(grid='True') jak tu zrobic boxplot???

print(rf2_best_params)
print("--------------------SVM (CV,GD)----------------------------------------")
clf_svm2 = svm.SVC(class_weight='balanced', random_state=22)
svm2_param_grid = {
    'kernel': ['linear', 'rbf', 'poly'],
    'C': [0.001, 0.01, 0.1, 1, 10],
    'gamma': [0.001, 0.01, 0.1, 1],
}
svm2_best_params = model_grid_search("svm", clf_svm2, x, y, svm2_param_grid)
clf_svm3 = svm.SVC(**svm2_best_params)
kfold(clf_svm3, x, y)
model_assessment("SVM", clf_svm3, x, y)
print(svm2_best_params)


def feature_plot(classifier, feature_names, top_features=20):
    coef = classifier.coef_.ravel()
    top_positive_coefficients = np.argsort(coef)[-top_features:]
    top_negative_coefficients = np.argsort(coef)[:top_features]
    top_coefficients = np.hstack([top_negative_coefficients, top_positive_coefficients])
    plt.figure(figsize=(18, 7))
    colors = ['green' if c < 0 else 'blue' for c in coef[top_coefficients]]
    plt.bar(np.arange(2 * top_features), coef[top_coefficients], color=colors)
    feature_names = np.array(feature_names)
    plt.xticks(np.arange(1 + 2 * top_features),
               feature_names[top_coefficients], rotation=45, ha='right')
    plt.show()


print(df.drop(['Outcome'], axis=1).columns.values)

trainedsvm = svm.LinearSVC().fit(X, Y)
feature_plot(trainedsvm, df.drop(['Outcome'], axis=1).columns.values)

print("--------------------XGBOOST (CV,GD)------------------------------------")
clf_xgb2 = XGBClassifier()
xgb2_param_grid = {
    'min_child_weight': [1, 4, 5, 6, 7, 9],
    'gamma': [0.5, 1, 1.5, 2, 5],
    'subsample': [0.5, 0.6, 0.7, 1.0],
    'colsample_bytree': [0.5, 0.6, 0.7, 1.0],
    'max_depth': [3, 4, 5, 6, 7, 10, 13, 15],
    'n_estimators': [100, 200, 300, 400, 500, 750, 1000]
}
xgboost_best_params = model_grid_search("XGBoost", clf_xgb2, x, y, xgb2_param_grid)
clf_xgb3 = XGBClassifier(**xgboost_best_params)
kfold(clf_xgb3, x, y)
model_assessment("XGBoost", clf_xgbc3, x, y)
print(xgboost_best_params)
