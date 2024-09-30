# import libraries
import os
import numpy as np
import matplotlib.pyplot as plt
import arviz as az
import pandas as pd
import hssm
plt.rcParams['figure.constrained_layout.use'] = True # overlap otherwise
hssm.set_floatX("float32")

# set cwd
os.chdir("C:\\Users\\ACER\\Desktop\\Codes\\DG\\data_anal")

# import data
import scipy.io
alldat = scipy.io.loadmat('behavDat_collated_Nov04.mat')
alldat = alldat['alldat']
print("data loaded")

# create dataset to fit ddm
dataset_ddm = []
k = 0;
for sub_idx in [12, 17]:
  dat = alldat[:,sub_idx][0]

  # skip for AA/AR
  accpt = sum(dat[:,-1]==1)/len(dat)
  if accpt > .95 or accpt < .05:
    print(sub_idx)

  else:
    # remove the missed trials
    dat = dat[dat[:,-1]!=0,:]
    # remove non-biological rt-s
    dat = dat[dat[:,3]>.3,:]

    for iter_idx in range(len(dat)):
      trial = dat[iter_idx, :]
      dataset_ddm.append([trial[3],
                          trial[-1],
                          trial[0]-1,
                          trial[1],
                          k])
    k += 1

dataset_ddm = pd.DataFrame(dataset_ddm,
                           columns = ['rt', 'response',
                                      'emot', 'off', 'subject'])
# scale offer
temp = dataset_ddm['off'];
dataset_ddm['off'] = (temp - temp.min()) / (temp.max() - temp.min())

# add dummy variable for emotion
dataset_ddm = pd.concat([dataset_ddm,
                         pd.get_dummies(dataset_ddm["emot"])],
                         axis = 1)

dataset_ddm = dataset_ddm.drop([0.0], axis = 1)
dataset_ddm = dataset_ddm.rename(columns = {'rt': 'rt',
                           'response': 'response',
                           'emot': 'emot_idx',
                           'off': 'off_norm',
                           'subject': 'subject',
                           1.0: 'emot_1',
                           2.0: 'emot_2'})
# keep the dummies as float
dataset_ddm[['emot_1', 'emot_2']] = dataset_ddm[['emot_1',
                                                 'emot_2']].astype(float)
dataset_ddm['emot_idx'] = dataset_ddm['emot_idx'].astype('category')

print("data table ready")

# define the ddm model
model_ddm_array = [];

for i in range(1):
  model_ddm_array.append(hssm.HSSM(
      data = dataset_ddm,
      model = 'ddm',
      include=[
          {
              "name": "v",
              "prior": {
                  "Intercept": {"name": "Normal", "mu": 0.0, "sigma": 2.0},
                  "off_norm": {"name": "Normal", "mu": 0.0, "sigma": 2.0},
                  "emot_idx": {"name": "Normal", "mu": 0.0, "sigma": 2.0},
              },
              "formula": "v ~ (1|subject) + off_norm + emot_idx",
              "link": "identity",
          },
          # {
          #     "name": "z",
          #     "prior": {
          #         "Intercept": {"name": "Uniform", "lower": 0, "upper": 1},
          #         "emot_idx": {"name": "Normal", "mu": 0.0, "sigma": 2.0},
          #     },
          #     "formula": "z ~ (1|subject) + emot_idx",
          #     "link": "logit",
          # },
          # {
          #     "name": "a",
          #     "prior": {
          #         "Intercept": {"name": "HalfNormal", "sigma": 2.0},
          #         "off_norm": {"name": "Normal", "mu": 0.0, "sigma": 2.0},
          #         "emot_idx": {"name": "Normal", "mu": 0.0, "sigma": 2.0},
          #     },
          #     "formula": "a ~ (1|subject) + off_norm + emot_idx",
          #     "link": "inverse_squared",
          # },
          # {
          #     "name": "t",
          #     "prior": {
          #         "Intercept": {"name": "HalfNormal", "sigma": 2.0},
          #         "off_norm": {"name": "Normal", "mu": 0.0, "sigma": 2.0},
          #         "emot_idx": {"name": "Normal", "mu": 0.0, "sigma": 2.0},
          #     },
          #     "formula": "t ~ (1|subject) + off_norm + emot_idx",
          #     "link": "inverse_squared",
          # }
      ],
  ))

print("model defined, now runs")
for i in range(1):
  idat_ddm = model_ddm_array[i].sample(chains=2,
                                      draws=1000,
                                      tune=1000,
                                      idata_kwargs=dict(
                                      log_likelihood=True),)
  
txt_file = open(r"C:\Users\ACER\Desktop\Codes\DG\data_anal\ddm_fit_results\res.txt")
txt_file.write(model_ddm_array[0].summary())
txt_file.close()

az.plot_trace(idat_ddm)
fig = plt.gcf()
fig.savefig(r"C:\Users\ACER\Desktop\Codes\DG\data_anal\ddm_fit_results\trace.png")

print("fig and result saved, done!")