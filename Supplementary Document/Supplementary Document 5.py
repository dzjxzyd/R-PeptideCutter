###### ABTS model
### input
import os
from sklearn.feature_selection import SelectFromModel
os.getcwd()
#### change your current working directory
os.chdir('/Users/zhenjiaodu/Desktop/QSAR_dipeptide')
import pandas as pd
import numpy as np
raw_ABTS= pd.read_excel('di_ABTS.xlsx') # read data

peptide_name = raw_ABTS['peptide_name'].to_numpy() # transform the peptide name into a numpy array
peptide_activty = raw_ABTS['activity'].to_numpy()# transform the peptide activity into a numpy array
#raw_ABTS
features= pd.read_csv('aa_553_washed_features.csv',header=0,index_col=0)
#features= pd.read_csv('aa_pearson_0.95_408_features.csv',header=0,index_col=0)
### transform peptide_name into a vector for model development
# creat a feature_for_model using the features to encode the peptide_name
sample_size = len(peptide_name) # get the peptide sampel size
peptide_length = len(peptide_name[0])
feature_dimension = np.shape(features)[0] # get dimension of the features
# create matirc for feature extraction
feature_for_model = np.empty([sample_size, peptide_length*feature_dimension])
np.shape(feature_for_model) # confirm the feature_for_model matrix dimenstion
for i in range(len(peptide_name)):
    name = peptide_name[i] # extract the peptide name; maybe a tripeptide
    try:
        first_aa = features[name[0]].to_numpy()
        second_aa = features[name[1]].to_numpy()
    except KeyError:
        pass
    # combine them together
    feature_for_model[i]= np.concatenate((first_aa,second_aa), axis=0)
raw_ABTS['activity'].shape
#### pearson screening
#transform the data type
feature_for_model = pd.DataFrame(feature_for_model)
Pearson_corr = feature_for_model.corr() # results is a pandas.core.frame.DataFrame
kk=Pearson_corr.to_numpy()
# make the Pearson_corr become a triangle matrix
for i in range(553):
    for j in range(i+1):
        kk[i][j]=0
kk=pd.DataFrame(kk)
kk.shape
#  get column name out
cols = list(Pearson_corr.columns)
selected_cols =  list(Pearson_corr.columns)
feature_for_model

for i in range(1106):
    for j in range(1106):
        if abs(Pearson_corr[cols[i]][j]) > 0.95: # set the coefficient for selection
            if i != j: # if not at the diagonal
                if cols[j] in selected_cols: # if still not be removed
                    selected_cols.remove(cols[j])
len(selected_cols)
# for the feature name obtaining;
# if beyond 553, then need to minus 553 first
selected_cols[0]
selected_cols[14]
selected_cols[99]
selected_cols[122]
selected_cols[130]
selected_cols[265]
selected_cols[376]
selected_cols[388]
selected_cols[466]
selected_cols[470]
selected_cols[521]
features.T.columns[selected_cols[0]]
features.T.columns[selected_cols[14]]
features.T.columns[selected_cols[99]]
features.T.columns[selected_cols[122]]
features.T.columns[selected_cols[130]]
features.T.columns[selected_cols[265]]
features.T.columns[selected_cols[376]]
features.T.columns[selected_cols[388]-553]
features.T.columns[selected_cols[466]-553]
features.T.columns[selected_cols[470]-553]
features.T.columns[selected_cols[521]-553]
features.T.shape
features.T['RICJ880105']

##array(['x0', 'x14', 'x99', 'x122', 'x130', 'x265', 'x376', 'x388', 'x466','x470', 'x521'], dtype=object)
#'original feature colomn is 22 168 393 472 672 687 878
len(selected_cols)
X = feature_for_model[selected_cols]
X.shape
y = peptide_activty
#df=pd.DataFrame(feature_for_model)
#df.to_csv('all_553*2.csv')
#df=pd.DataFrame(X)
#df.to_csv('after selected.csv')
mean = X.mean().to_numpy()
std = X.std().to_numpy()
#(X-X.mean()).to_numpy() / X.std().to_numpy()

# standardlization
from sklearn.preprocessing import StandardScaler
sc = StandardScaler() #
X= sc.fit_transform(X)
from dataset_split import dataset_split
X_train,X_test, y_train,y_test = dataset_split(X,y)
# XGBRegressor +
import pandas as pd
# XGB for selection
import xgboost as xgb
rf = xgb.XGBRegressor(max_depth=2,n_estimators=800,reg_alpha= 2,  reg_lambda = 0) # set the l2 = 1 to avoid over fitting
rf.fit(X_train, y_train)
rf.score(X_test,y_test)
importances = rf.feature_importances_ # extract the feature_importances_

# np.argsort(importances)
sfm = SelectFromModel(rf, threshold=0.01)
sfm.fit(X_train, y_train)
# get the selected features and thier feature_importance values
sfm.get_feature_names_out()

rf.feature_importances_[0]
rf.feature_importances_[14]
rf.feature_importances_[99]
rf.feature_importances_[122]
rf.feature_importances_[130]
rf.feature_importances_[265]
rf.feature_importances_[376]
rf.feature_importances_[388]
rf.feature_importances_[466]
rf.feature_importances_[470]
rf.feature_importances_[521]
#array(['x0', 'x14', 'x99', 'x122', 'x130', 'x265', 'x376', 'x388', 'x466','x470', 'x521'], dtype=object)
#note: the X12 means the 13th features,

#rf = RandomForestRegressor(n_estimators = 1000,  random_state=0)
#rf.fit(X_important_train, y_train)
#rf.score(X_important_train, y_train)
#rf.score(X_important_test,y_test)

#####
#####
##### data transform and shuffle
X_new = sfm.transform(X)
X_new.shape
# dateset seperation, randomly division
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle
X_shuffle, y_shuffle = shuffle(X_new, y,random_state=0)
X_important_train, X_important_test, y_train, y_test = train_test_split(X_shuffle, y_shuffle, test_size=0.25,random_state = 1)
#### test begin
#### 1. xgbtree
import xgboost as xgb
from sklearn.metrics import r2_score, mean_squared_error, explained_variance_score
##### cross validation
rf = xgb.XGBRegressor(max_depth=2,n_estimators=800,reg_alpha= 2.5,  reg_lambda = 0)
#(max_depth=2,n_estimators=800,reg_alpha= 3,  reg_lambda = 1.1)  for 8 feastures
#
model = rf
model.fit(X_important_train, y_train)
train_score = model.score(X_important_train, y_train)
train_score
y_train_pre = model.predict(X_important_train)
train_MSE= mean_squared_error(y_train, y_train_pre)
train_MSE
test_score= model.score(X_important_test, y_test)
test_score
y_test_pre = model.predict(X_important_test)
test_MSE=mean_squared_error(y_test, y_test_pre)
test_MSE
y_cross_val_pred = []
y_cross_val_test=[]
for i in range(X_important_train.shape[0]):
    rest_X_train = np.delete(X_important_train,i,axis=0)
    rest_y_train = np.delete(y_train,i,axis=0)
    out_X = np.empty(([1,11]))
    out_X[0][0] = X_important_train[i][0]
    out_X[0][1] = X_important_train[i][1]
    out_X[0][2] = X_important_train[i][2]
    out_X[0][3] = X_important_train[i][3]
    out_X[0][4] = X_important_train[i][4]
    out_X[0][5] = X_important_train[i][5]
    out_X[0][6] = X_important_train[i][6]
    out_X[0][7] = X_important_train[i][7]
    out_X[0][8] = X_important_train[i][8]
    out_X[0][9] = X_important_train[i][9]
    #out_X[0][10] = X_important_train[i][10]
    #out_X[0][11] = X_important_train[i][11]
    model.fit(rest_X_train, rest_y_train)
    results = model.predict(out_X)
    y_cross_val_pred.append(results[0])
    y_cross_val_test.append(y_train[i])
from sklearn.metrics import r2_score, mean_squared_error
r2_cross_val=r2_score(y_cross_val_test, y_cross_val_pred)
r2_cross_val
mse_cross_val=mean_squared_error(y_cross_val_test, y_cross_val_pred)
mse_cross_val






####prediction
####
pred_ABTS= pd.read_excel('all the possible diprptide.xlsx') # read data
pred_name = pred_ABTS['peptide_name'].to_numpy() # transform the peptide name into a numpy array
features= pd.read_csv('aa_553_washed_features.csv',header=0,index_col=0)
sample_size = len(pred_name) # get the peptide sampel size
peptide_length = len(pred_name[0])
feature_dimension = np.shape(features)[0] # get dimension of the features
# create matirc for feature extraction
feature_for_model = np.empty([sample_size, peptide_length*feature_dimension])
np.shape(feature_for_model) # confirm the feature_for_model matrix dimenstion
for i in range(len(pred_name)):
    name = pred_name[i] # extract the peptide name; maybe a tripeptide
    if i != 40:
        try:
            first_aa = features[name[0]].to_numpy()
            second_aa = features[name[1]].to_numpy()
        except KeyError:
            pass
    if i == 40:
        name = 'NA'
        try:
            first_aa = features[name[0]].to_numpy()
            second_aa = features[name[1]].to_numpy()
        except KeyError:
            pass
    # combine them together
    feature_for_model[i]= np.concatenate((first_aa,second_aa), axis=0)
feature_for_model.shape

# used the selected feature by Pearson to find the feature used in the unknown data
pred_X = feature_for_model.T[selected_cols]
pred_X.shape
pred_X= pred_X.T
pred_X.shape
# standardlization bsed on the method we used in previous dataset
pred_X = (pred_X-mean) / std #
pre_feature = sfm.transform(pred_X)
pre_feature.shape
pre_activity = model.predict(pre_feature)

df = pd.DataFrame(pre_activity)
df.to_csv('ABTS_potential activity results_new_for_scaler1.csv')

# find the published values
for i in range(len(pred_name)):
    name = pred_name[i] #
    if name in peptide_name:
        pre_activity[i] = 6666 # set a strange value to distinguish the published values
# sorted the results as a matrix
sorted_pep= dict()
for i in range(20):
    sorted_pep[i]=pre_activity[20*i:(i+1)*20]
df_sorted = pd.DataFrame(sorted_pep) # column name is the N_amino acid
df_sorted.to_csv('matrix_ABTS_potential activity results.csv')












#####  ORAC model
### input
import os
from sklearn.feature_selection import SelectFromModel
os.getcwd()
#### change your current working directory
os.chdir('/Users/zhenjiaodu/Desktop/QSAR_dipeptide')
import pandas as pd
import numpy as np
raw_ORAC= pd.read_excel('di_ORAC.xlsx') # read data

peptide_name = raw_ORAC['peptide_name'].to_numpy() # transform the peptide name into a numpy array
peptide_activty = raw_ORAC['activity'].to_numpy()# transform the peptide activity into a numpy array
#raw_ORAC
features= pd.read_csv('aa_553_washed_features.csv',header=0,index_col=0)
#features= pd.read_csv('aa_pearson_0.95_408_features.csv',header=0,index_col=0)
### transform peptide_name into a vector for model development
# creat a feature_for_model using the features to encode the peptide_name
sample_size = len(peptide_name) # get the peptide sampel size
peptide_length = len(peptide_name[0])
feature_dimension = np.shape(features)[0] # get dimension of the features
# create matirc for feature extraction
feature_for_model = np.empty([sample_size, peptide_length*feature_dimension])
np.shape(feature_for_model) # confirm the feature_for_model matrix dimenstion
for i in range(len(peptide_name)):
    name = peptide_name[i] # extract the peptide name; maybe a tripeptide
    try:
        first_aa = features[name[0]].to_numpy()
        second_aa = features[name[1]].to_numpy()
    except KeyError:
        pass
    # combine them together
    feature_for_model[i]= np.concatenate((first_aa,second_aa), axis=0)
raw_ORAC['activity'].shape
#### pearson screening
#transform the data type
feature_for_model = pd.DataFrame(feature_for_model)
Pearson_corr = feature_for_model.corr() # results is a pandas.core.frame.DataFrame
kk=Pearson_corr.to_numpy()
# make the Pearson_corr become a triangle matrix
for i in range(553):
    for j in range(i+1):
        kk[i][j]=0
kk=pd.DataFrame(kk)
kk.shape
#  get column name out
cols = list(Pearson_corr.columns)
selected_cols =  list(Pearson_corr.columns)
# get the feature name that beyond Pearson coefficient above 0.95
for i in range(1106):
    for j in range(1106):
        if abs(Pearson_corr[cols[i]][j]) > 0.95: # set the coefficient for selection
            if i != j: # if not at the diagonal
                if cols[j] in selected_cols: # if still not be removed
                    selected_cols.remove(cols[j])
len(selected_cols)
selected_cols[5]
selected_cols[6]
selected_cols[23]
selected_cols[32]
selected_cols[100]
selected_cols[131]
selected_cols[310]
selected_cols[388]
selected_cols[396]
selected_cols[407]
selected_cols[417]
selected_cols[498]
selected_cols[539]
selected_cols[674]

features.T.columns[selected_cols[5]]
features.T.columns[selected_cols[6]]
features.T.columns[selected_cols[23]]
features.T.columns[selected_cols[32]]
features.T.columns[selected_cols[100]]
features.T.columns[selected_cols[131]]
features.T.columns[selected_cols[310]]
features.T.columns[selected_cols[388]-553]
features.T.columns[selected_cols[396]-553]
features.T.columns[selected_cols[407]-553]
features.T.columns[selected_cols[417]-553]
features.T.columns[selected_cols[498]-553]
features.T.columns[selected_cols[539]-553]
features.T.columns[selected_cols[674]-553]
features.T.shape
#array(['x5', 'x6', 'x23', 'x32', 'x100', 'x131', 'x310', 'x388', 'x396',
#       'x407', 'x417', 'x498', 'x539', 'x674'], dtype=object)
#'original feature colomn is 22 168 393 472 672 687 878
len(selected_cols)
X = feature_for_model[selected_cols]
X.shape
y = peptide_activty
#df=pd.DataFrame(feature_for_model)
#df.to_csv('all_553*2.csv')
#df=pd.DataFrame(X)
#df.to_csv('after selected.csv')
mean = X.mean().to_numpy()
std = X.std().to_numpy()
#(X-X.mean()).to_numpy() / X.std().to_numpy()

# standardlization
from sklearn.preprocessing import StandardScaler
sc = StandardScaler() #
X= sc.fit_transform(X)
X.shape
from dataset_split import dataset_split
X_train,X_test, y_train,y_test = dataset_split(X,y)
# XGBRegressor +
import pandas as pd
# XGB for selection
import xgboost as xgb
rf = xgb.XGBRegressor(max_depth=2,n_estimators=800,reg_alpha= 2.8,  reg_lambda = 8.1) # set the l2 = 1 to avoid over fitting
rf.fit(X_train, y_train)
rf.score(X_train, y_train)
rf.score(X_test,y_test)

importances = rf.feature_importances_ # extract the feature_importances_

# np.argsort(importances)
sfm = SelectFromModel(rf, threshold=0.01)
sfm.fit(X_train, y_train)
# get the selected features and thier feature_importance values
sfm.get_feature_names_out()

#note: the X12 means the 13th features,
rf.feature_importances_[5]
rf.feature_importances_[6]
rf.feature_importances_[23]
rf.feature_importances_[32]
rf.feature_importances_[100]
rf.feature_importances_[131]
rf.feature_importances_[310]
rf.feature_importances_[388]
rf.feature_importances_[396]
rf.feature_importances_[407]
rf.feature_importances_[417]
rf.feature_importances_[498]
rf.feature_importances_[539]
rf.feature_importances_[674]


#####
##### data transform and shuffle
X_new = sfm.transform(X)
X_new.shape
# dateset seperation, randomly division
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle
X_shuffle, y_shuffle = shuffle(X_new, y,random_state=0)
X_important_train, X_important_test, y_train, y_test = train_test_split(X_shuffle, y_shuffle, test_size=0.25,random_state = 0)
#### test begin
#### 1. xgbtree
import xgboost as xgb
from sklearn.metrics import r2_score, mean_squared_error, explained_variance_score
##### cross validation
rf = xgb.XGBRegressor(max_depth=1,n_estimators=500,reg_alpha= 0.9,  reg_lambda =1.8)
model = rf
model.fit(X_important_train, y_train)
train_score = model.score(X_important_train, y_train)
train_score
y_train_pre = model.predict(X_important_train)
train_MSE= mean_squared_error(y_train, y_train_pre)
train_MSE
test_score= model.score(X_important_test, y_test)
test_score
y_test_pre = model.predict(X_important_test)
test_MSE=mean_squared_error(y_test, y_test_pre)
test_MSE
y_cross_val_pred = []
y_cross_val_test=[]
for i in range(X_important_train.shape[0]):
    rest_X_train = np.delete(X_important_train,i,axis=0)
    rest_y_train = np.delete(y_train,i,axis=0)
    out_X = np.empty(([1,14]))
    out_X[0][0] = X_important_train[i][0]
    out_X[0][1] = X_important_train[i][1]
    out_X[0][2] = X_important_train[i][2]
    out_X[0][3] = X_important_train[i][3]
    out_X[0][4] = X_important_train[i][4]
    out_X[0][5] = X_important_train[i][5]
    out_X[0][6] = X_important_train[i][6]
    out_X[0][7] = X_important_train[i][7]
    out_X[0][8] = X_important_train[i][8]
    out_X[0][9] = X_important_train[i][9]
    out_X[0][10] = X_important_train[i][10]
    out_X[0][11] = X_important_train[i][11]
    out_X[0][12] = X_important_train[i][12]
    out_X[0][13] = X_important_train[i][13]
    model.fit(rest_X_train, rest_y_train)
    results = model.predict(out_X)
    y_cross_val_pred.append(results[0])
    y_cross_val_test.append(y_train[i])
from sklearn.metrics import r2_score, mean_squared_error
r2_cross_val=r2_score(y_cross_val_test, y_cross_val_pred)
r2_cross_val
mse_cross_val=mean_squared_error(y_cross_val_test, y_cross_val_pred)
mse_cross_val


#### prediction
pred_ORAC= pd.read_excel('all the possible diprptide.xlsx') # read data
pred_name = pred_ORAC['peptide_name'].to_numpy() # transform the peptide name into a numpy array
features= pd.read_csv('aa_553_washed_features.csv',header=0,index_col=0)
sample_size = len(pred_name) # get the peptide sampel size
peptide_length = len(pred_name[0])
feature_dimension = np.shape(features)[0] # get dimension of the features
# create matirc for feature extraction
feature_for_model = np.empty([sample_size, peptide_length*feature_dimension])
np.shape(feature_for_model) # confirm the feature_for_model matrix dimenstion
for i in range(len(pred_name)):
    name = pred_name[i] # extract the peptide name; maybe a tripeptide
    if i != 40:
        try:
            first_aa = features[name[0]].to_numpy()
            second_aa = features[name[1]].to_numpy()
        except KeyError:
            pass
    if i == 40:
        name = 'NA'
        try:
            first_aa = features[name[0]].to_numpy()
            second_aa = features[name[1]].to_numpy()
        except KeyError:
            pass
    # combine them together
    feature_for_model[i]= np.concatenate((first_aa,second_aa), axis=0)
feature_for_model.shape

# used the selected feature by Pearson to find the feature used in the unknown data
pred_X = feature_for_model.T[selected_cols]
pred_X.shape
pred_X= pred_X.T
pred_X.shape
# standardlization bsed on the method we used in previous dataset
pred_X = (pred_X-mean) / std #
pre_feature = sfm.transform(pred_X)
pre_activity = model.predict(pre_feature)

df = pd.DataFrame(pre_activity)
df.to_csv('ORAC_potential activity results_new_for_scaler1.csv')

# find the published values
for i in range(len(pred_name)):
    name = pred_name[i] #
    if name in peptide_name:
        pre_activity[i] = 6666 # set a strange value to distinguish the published values
# sorted the results as a matrix
sorted_pep= dict()
for i in range(20):
    sorted_pep[i]=pre_activity[20*i:(i+1)*20]
df_sorted = pd.DataFrame(sorted_pep) # column name is the N_amino acid
df_sorted.to_csv('matrix_ORAC_potential activity results as 6666.csv')
