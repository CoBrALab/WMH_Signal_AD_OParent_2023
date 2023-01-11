#%% Initialization

from pyls import behavioral_pls
import numpy as np
import pandas as pd
import os
import sys
from numpy import genfromtxt

# Set working directory to script location
wd = os.path.dirname(os.path.realpath(__file__))
os.chdir(wd)

#%% Make X and Y matrices

ADB_subset = pd.read_csv("../ADB_subset_imputed_RF_new.csv")
ADB_subset = ADB_subset[ADB_subset.columns.drop(list(ADB_subset.filter(regex='NAWM')))]
ADB_subset = ADB_subset.iloc[:,:-4]
ADB_subset['TBV_ICV_ratio'] = ADB_subset['TBV_ICV_ratio'] / 1000
ADB_subset['TBV'] = ADB_subset['TBV'] / 1000

list(ADB_subset.columns)

# Cognition: RBANS subscales and total
Y_raw = ADB_subset[['Age', 'Sex', 'BMI', 'Years_school', 'MoCA_total', 'AD8', 'APOE4_status', 
                    "RBANS_immediate_memory", "RBANS_visuospatial_memory", "RBANS_language", "RBANS_attention", "RBANS_delayed_memory", "RBANS_total",
                    "High_BP", "High_chol", "Alcohol"]]
# Brain: Global WMH characteristics (no parcellation)
X_raw = ADB_subset.iloc[:,20:]
#X_raw = X_raw.loc[:,~X_raw.columns.isin(['TBV', 'ICV'])]

# Check for missing data
print(Y_raw.isnull().sum(axis = 0))
print(X_raw.isnull().sum(axis = 0))

y_nan = Y_raw.isnull().sum(axis = 0)
x_nan = X_raw.isnull().sum(axis = 0)

# Check variance
print(Y_raw.var())
print(X_raw.var())

# Convert to numpy array
Y = Y_raw.to_numpy()
X = X_raw.to_numpy()

print(np.isnan(np.sum(Y)))
print(np.isnan(np.sum(X)))

# Create a directory in which to save the outputs
output_dir = (f'{wd}/outputs/')
os.makedirs(output_dir, exist_ok = True)

#%% Run PLS and save outputs

# Run plsc
bpls = behavioral_pls(X,Y, n_perm = 5000, n_boot = 5000, n_split = 200, seed = 123)
# n_proc can be set to ‘max’ on niagara

np.savetxt(f'{output_dir}/x_weights.csv', bpls['x_weights'], delimiter=',') # p x m
np.savetxt(f'{output_dir}/y_weights.csv', bpls['y_weights'], delimiter=',') # m x m
np.savetxt(f'{output_dir}/x_scores.csv', bpls['x_scores'], delimiter=',') # n x m
np.savetxt(f'{output_dir}/y_scores.csv', bpls['y_scores'], delimiter=',') # n x m
np.savetxt(f'{output_dir}/y_loadings.csv', bpls['y_loadings'], delimiter=',') # m x m 
np.savetxt(f'{output_dir}/sinvals.csv', bpls['singvals'], delimiter=',') # (m,)
np.savetxt(f'{output_dir}/varexp.csv', bpls['varexp'], delimiter=',') # (m,)
np.savetxt(f'{output_dir}/permres_pvals.csv', bpls['permres']['pvals'], delimiter=',') # shape: (m,)
np.savetxt(f'{output_dir}/permres_permsamples.csv', bpls['permres']['permsamples'], delimiter=',') 
np.savetxt(f'{output_dir}/bootres_x_weights_normed.csv', bpls['bootres']['x_weights_normed'], delimiter=',') 
np.savetxt(f'{output_dir}/bootres_x_weights_stderr.csv', bpls['bootres']['x_weights_stderr'], delimiter=',')
np.savetxt(f'{output_dir}/bootres_bootsamples.csv', bpls['bootres']['bootsamples'], delimiter=',')
np.savetxt(f'{output_dir}/bootres_y_loadings.csv', bpls['bootres']['y_loadings'], delimiter=',')

# to save your bootstrapped samples (for the very, very cautious)
os.makedirs(f'{output_dir}/y_loadings_boot', exist_ok = True)
for i in range(bpls['bootres']['y_loadings_boot'].shape[0]):
    np.savetxt(f'{output_dir}/y_loadings_boot/bootres_y_loadings_boot_{i}.csv', bpls['bootres']['y_loadings_boot'][i], delimiter=',')

# saving your confidence intervals 
os.makedirs(f'{output_dir}/y_loadings_ci', exist_ok = True)
for i in range(bpls['bootres']['y_loadings_ci'].shape[0]):
    np.savetxt(f'{output_dir}/y_loadings_ci/bootres_y_loadings_ci_behaviour_{i}.csv', bpls['bootres']['y_loadings_ci'][i], delimiter=',')

np.savetxt(f'{output_dir}/cvres_pearson_r.csv', bpls['cvres']['pearson_r'], delimiter=',') 
np.savetxt(f'{output_dir}/cvres_r_squared.csv', bpls['cvres']['r_squared'], delimiter=',')

# Saving your input information (makes the code longer, but quite valuable when trouble shooting, or if you've done several runs and want to keep track of the parameters of each):
np.savetxt(f'{output_dir}/inputs_X.csv', bpls['inputs']['X'], delimiter=',')
np.savetxt(f'{output_dir}/inputs_Y.csv', bpls['inputs']['Y'], delimiter=',')
f=open(f'{output_dir}/input_info.txt','w')
f.write(f"Groups: {bpls['inputs']['groups']}\n") # [n]
f.write(f"Num of conditions: {bpls['inputs']['n_cond']}\n") # 1
f.write(f"Num of permutations: {bpls['inputs']['n_perm']}\n") # perm
f.write(f"Bootstrapping: {bpls['inputs']['n_boot']}\n") # boot
f.write(f"Bootstrapping: {bpls['inputs']['test_split']}\n") # default=100
f.write(f"Bootstrapping: {bpls['inputs']['test_size']}\n") # default=0.25
f.write(f"Bootstrapping: {bpls['inputs']['covariance']}\n") # default=False
f.write(f"Rotations: {bpls['inputs']['rotate']}\n") # True
f.write(f"Confidence Intervals: {bpls['inputs']['ci']}\n") # default=95
f.write(f"Verbose: {bpls['inputs']['verbose']}\n") # True
f.close()

# Saving split-half resampling results
os.makedirs(f'{output_dir}/splitres/', exist_ok=True)
np.savetxt(f'{output_dir}/splitres/splitres_u-vcorr.csv', np.column_stack((bpls['splitres']['ucorr'], bpls['splitres']['vcorr'])), delimiter=',')
np.savetxt(f'{output_dir}/splitres/splitres_u-vcorr_pvals.csv', np.column_stack((bpls['splitres']['ucorr_pvals'], bpls['splitres']['vcorr_pvals'])), delimiter=',')
np.savetxt(f'{output_dir}/splitres/splitres_ucorr_lo-uplim.csv', np.column_stack((bpls['splitres']['ucorr_lolim'], bpls['splitres']['ucorr_uplim'])),delimiter=',')
np.savetxt(f'{output_dir}/splitres/splitres_vcorr_lo-uplim.csv', np.column_stack((bpls['splitres']['vcorr_lolim'], bpls['splitres']['vcorr_uplim'])),delimiter=',')


#%% Import results (instead of rerunning)
import shutil

results_dir=output_dir # directory in which your PLS outputs have been saved
shutil.rmtree(f"{wd}/visualization", ignore_errors=True)
os.makedirs(f"{wd}/visualization", exist_ok = True)
vis_dir=(f"{wd}/visualization") # directory in which you will be saving your PLS visualizations
df = Y_raw # data frame where each row is a subject, each column is a behavioural variable

inputs_X = pd.read_csv(f"{results_dir}/inputs_X.csv", header=None)
inputs_Y = pd.read_csv(f"{results_dir}/inputs_Y.csv", header=None)
print("inputs_X shape: ",inputs_X.shape)
print("inputs_Y shape: ",inputs_Y.shape)
x_scores = pd.read_csv(f"{results_dir}/x_scores.csv", header=None)
y_scores = pd.read_csv(f"{results_dir}/y_scores.csv", header=None)
print("x_scores shape: ",x_scores.shape)
print("y_scores shape: ",y_scores.shape)
bootres_bootsamples = pd.read_csv(f"{results_dir}/bootres_bootsamples.csv", header=None)
print("bootres_bootsamples shape: ",bootres_bootsamples.shape)
bootres_x_weights_normed = pd.read_csv(f"{results_dir}/bootres_x_weights_normed.csv", header=None)
print("bootres_x_weights_normed shape: ",bootres_x_weights_normed.shape)
bootres_x_weights_stderr = pd.read_csv(f"{results_dir}/bootres_x_weights_stderr.csv", header=None)
print("bootres_x_weights_stderr shape: ",bootres_x_weights_stderr.shape)
permres_permsamples = pd.read_csv(f"{results_dir}/permres_permsamples.csv", header=None)
print("permres_permsamples shape: ",permres_permsamples.shape)
permres_pvals= pd.read_csv(f"{results_dir}/permres_pvals.csv", header=None)
print("permres_pvals shape: ",permres_pvals.shape)
varexp = pd.read_csv(f"{results_dir}/varexp.csv", header=None)
print("varexp shape: ",varexp.shape)
x_weights= pd.read_csv(f"{results_dir}/x_weights.csv", header=None)
print("x_weights shape: ",x_weights.shape)
y_loadings= pd.read_csv(f"{results_dir}/y_loadings.csv", header=None)
print("y_loadings shape: ",y_loadings.shape)
splitres_u_vcorr_pvals = pd.read_csv(f"{results_dir}splitres/splitres_u-vcorr_pvals.csv", header=None)
print("splitres_u_vcorr_pvals shape: ",splitres_u_vcorr_pvals.shape)


#%% Visualization: P-values and covariance explained by each latent variable
import math
import matplotlib.pyplot as plt

# P-values and covariance explained by each latent variable
fig, ax1 = plt.subplots()
LV = list(range(1,df.shape[1]+1))

color = 'tab:red'
ax1.set_xlabel('LV')

ax1.set_ylabel('% covariance', color=color)
ax1.scatter(LV, varexp*100, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.tick_params(axis='x')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('p-value', color=color)  # we already handled the x-label with ax1
ax2.scatter(LV, permres_pvals, color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.axhline(y=0.05, color="tab:blue")
ax1.xaxis.grid(which='major')

fig.tight_layout()  # otherwise the right y-label is slightly clipped

ax1.set_title("Covariance explained and p-values")

ax1.title.set_fontsize(20)
ax1.xaxis.label.set_fontsize(20)
ax1.yaxis.label.set_fontsize(20)
ax2.yaxis.label.set_fontsize(20)
ax1.tick_params(axis='x', labelsize=15)
ax1.tick_params(axis='y', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])

#fig.savefig(f"{vis_dir}/pval_by_covexp_per_LV.png")
fig.savefig(f"{vis_dir}/pval_by_covexp_per_LV.png",bbox_inches = "tight", dpi=1000)

#%% Visualization: Split-half resampling p-values
import math
import matplotlib.pyplot as plt

# P-values and covariance explained by each latent variable
fig, ax1 = plt.subplots()
LV = list(range(1,df.shape[1]+1))

LV_sp1=list()
LV_sp2=list()
LV_sp1[:] = [i-0.1 for i in LV]
LV_sp2[:] = [i+0.1 for i in LV]

color = 'tab:red'
ax1.set_xlabel('LV')

ax1.set_ylabel('P Ucorr', color=color)
ax1.scatter(LV_sp1, np.minimum(splitres_u_vcorr_pvals[0], 0.199), color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.tick_params()
ax1.set_ylim([0, 0.20])

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('P Vcorr', color=color)  # we already handled the x-label with ax1
ax2.scatter(LV_sp2, np.minimum(splitres_u_vcorr_pvals[1], 0.199), color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim([0, 0.20])

plt.axhline(y=0.05, color="tab:blue")
ax1.xaxis.grid(which='major')

fig.tight_layout()  # otherwise the right y-label is slightly clipped

ax1.set_title("Split-half resampling p-values")

ax1.title.set_fontsize(20)
ax1.xaxis.label.set_fontsize(20)
ax1.yaxis.label.set_fontsize(20)
ax2.yaxis.label.set_fontsize(20)
ax1.tick_params(axis='x', labelsize=15)
ax1.tick_params(axis='y', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])

fig.savefig(f"{vis_dir}/splitres_u_vcorr_pvals.png",bbox_inches = "tight", dpi=1000)


#%% Visualization: behavioral loadings on LVs
behav = list(Y_raw.columns)
y_pos = list(range(0,len(behav)))

for i in range(permres_pvals.shape[0]): # i ranges from 0 to (num_LVs – 1)
    if (permres_pvals[0][i] <= 0.05):  
        y_load = y_loadings.iloc[i,:]
        x_error = []
        for j in range(len(behav)):
            y_loadings_ci = np.loadtxt(f"{results_dir}y_loadings_ci/bootres_y_loadings_ci_behaviour_{j}.csv",delimiter=",")
            pair = y_loadings_ci[i]
            entry = [math.fabs(float((y_load[j]-pair[0]))),math.fabs(float((y_load[j]-pair[1])))]
            x_error.append(entry)
            
        # SORTING LOADINGS IN ASCENDING ORDER BY:
        # x_error_new = [x_error[i] for i in new_indices]
        # colors_new = [colors[i] for i in new_indices]
        # behav_new = [behav[i] for i in new_indices]
        # y_load.sort()
        # x_error_new = np.array(x_error_new)
        # x_error_new = x_error_new.T
        x_error = np.array(x_error)
        x_error = x_error.T
        
        # COLOURING IT SUCH THAT ONLY THE BEHAVIOURS WITH LOWER & UPPER CIs IN THE SAME DIRECTION APPEAR IN RED:
        colors=[]
        for k in range(x_error.shape[1]):
            a=y_load[k]-x_error[0,k]
            b=y_load[k]+x_error[1,k]
            c=a*b
            if (c >= 0):
                colors.append("green")
            elif (c < 0):
                colors.append("grey")
        
        plt.rc('font', size=30)
        fig, ax = plt.subplots(figsize=(15, 10),dpi=300)
        y_load = y_load*-1 # Invert loadings to facilitate interpretation
        x_error = x_error*-1 # Invert loadings to facilitate interpretation
        ax.barh(y_pos, y_load,xerr=x_error, align='center',color=colors)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(behav)
        for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), colors):
            ticklabel.set_color(tickcolor)
        ax.invert_yaxis()
        ax.set_xlabel('Loadings')
        ax.set_title(f"LV{i+1} Behavioural Loadings")
        ax.yaxis.grid()
        #plt.tight_layout()
        plt.savefig(f"{vis_dir}/LV_{i+1}_behav", bbox_inches='tight')
        plt.show()    

#%% Visualization: Brain loadings on LVs
for i in range(permres_pvals.shape[0]): # i ranges from 0 to (num_LVs – 1)
    if (permres_pvals[0][i] <= 0.05):  
        brain = list(X_raw.columns)
        y_pos = list(range(0,len(brain)))
        y_load = bootres_x_weights_normed.iloc[:,i]
        y_load = y_load*-1 # Invert BSR to facilitate interpretation
        y_load_abs = abs(bootres_x_weights_normed.iloc[:,i])

        # Sort by BSR
        y_all = pd.DataFrame([brain, y_pos, y_load, y_load_abs]).T
        y_all.columns = ['brain', 'y_pos', 'y_load', 'y_load_abs']
        y_all = y_all.sort_values(by=['y_load_abs'], ascending=False).reset_index()
    
        BSR_thresh_med = 2.57 #p=0.001;BSR=3.29  p=0.01;BSR=2.57  p=0.05;BSR=1.95
        BSR_thresh_high= 3.29 #p=0.001;BSR=3.29  p=0.01;BSR=2.57  p=0.05;BSR=1.95

        # Color: red if above BSR threshold, blue otherwise
        colors=[]
        for k in range(len(y_all['y_load'])):
            c=y_all['y_load'][k]
            if (abs(c) >= BSR_thresh_high):
                if ('CT' in y_all['brain'][k]):
                    colors.append("red")
                elif ('TBV' in y_all['brain'][k]):
                    colors.append("red")
                elif ('ICV' in y_all['brain'][k]):
                    colors.append("red")
                else:
                    colors.append("blue")

            #elif (abs(c) >= BSR_thresh_med and abs(c) < BSR_thresh_high):
            #    colors.append("red")
            elif (abs(c) < BSR_thresh_high):
                colors.append("grey")
        
        plt.rc('font', size=10)
        fig, ax = plt.subplots(figsize=(10, 10),dpi=300)
        ax.barh(y_all['y_load'].index, y_all['y_load'], align='center', color=colors)
        ax.set_yticks(y_all['y_load'].index)
        ax.set_yticklabels(y_all['brain'])
        for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), colors):
            ticklabel.set_color(tickcolor)
        ax.invert_yaxis()
        ax.set_xlabel('Bootstrap ratio')
        ax.set_title(f"LV{i+1} Brain Loadings")
        #plt.axvline(x=BSR_thresh_med, color="red")
        #plt.axvline(x=-BSR_thresh_med, color="red")
        plt.axvline(x=BSR_thresh_high, color="black")
        plt.axvline(x=-BSR_thresh_high, color="black")
        ax.yaxis.grid()
        #plt.tight_layout()
        plt.savefig(f"{vis_dir}/LV_{i+1}_brain", bbox_inches='tight')
        plt.show()    