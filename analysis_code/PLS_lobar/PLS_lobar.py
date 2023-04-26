# Partial Least Squares Correlation (PLSC) analysis relating multivariate
# patterns of brain variables (WMH measures and atrophy) to cognitive,
# demographic, and cardiovascular risk factor variables

# Implemented with the pyls package (https://pyls.readthedocs.io/en/latest/generated/pyls.behavioral_pls.html)

# Parcellation: Lobar

#%% Initialization

from pyls import behavioral_pls
import numpy as np
import pandas as pd
import os
from numpy import genfromtxt
import math
import matplotlib.pyplot as plt

# Set working directory to script location
wd = os.path.dirname(os.path.realpath(__file__))
os.chdir(wd)

#%% Make X and Y matrices

ADB_subset = pd.read_csv("../../df_lobar.csv")
ADB_subset = ADB_subset[ADB_subset.columns.drop(list(ADB_subset.filter(regex='NAWM')))]
ADB_subset['TBV_ICV_ratio'] = ADB_subset['TBV_ICV_ratio'] / 1000
ADB_subset['TBV'] = ADB_subset['TBV'] / 1000

list(ADB_subset.columns)

# Y matrix: Demographics, cognition, cardiovascular risk factors
Y_raw = ADB_subset[['Age', 'Sex', 'BMI', 'Years_school', 'MoCA_total', 'AD8', 'APOE4_status',
                    "RBANS_immediate_memory", "RBANS_visuospatial_memory", "RBANS_language", "RBANS_attention", "RBANS_delayed_memory", "RBANS_total",
                    "High_BP", "High_chol", "Alcohol"]]

# X matrix: WMH and atrophy measures
X_raw = ADB_subset.iloc[:,20:]

# Convert to numpy array
Y = Y_raw.to_numpy()
X = X_raw.to_numpy()

# Create a directory in which to save the outputs
output_dir = (f'{wd}/results/')
os.makedirs(output_dir, exist_ok = True)

#%% Run PLS and save outputs

# Run plsc
bpls = behavioral_pls(X,Y, n_perm = 5000, n_boot = 5000, n_split = 200, seed = 527)

# Save results
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

# Saving your confidence intervals 
os.makedirs(f'{output_dir}/y_loadings_ci', exist_ok = True)
for i in range(bpls['bootres']['y_loadings_ci'].shape[0]):
    np.savetxt(f'{output_dir}/y_loadings_ci/bootres_y_loadings_ci_behaviour_{i}.csv', bpls['bootres']['y_loadings_ci'][i], delimiter=',')

np.savetxt(f'{output_dir}/cvres_pearson_r.csv', bpls['cvres']['pearson_r'], delimiter=',') 
np.savetxt(f'{output_dir}/cvres_r_squared.csv', bpls['cvres']['r_squared'], delimiter=',')

# Saving your input information
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

#%% Import results

results_dir=output_dir
os.makedirs(f"{wd}/visualization", exist_ok = True)
vis_dir=(f"{wd}/visualization")
df = Y_raw

inputs_X = pd.read_csv(f"{results_dir}/inputs_X.csv", header=None)
inputs_Y = pd.read_csv(f"{results_dir}/inputs_Y.csv", header=None)
x_scores = pd.read_csv(f"{results_dir}/x_scores.csv", header=None)
y_scores = pd.read_csv(f"{results_dir}/y_scores.csv", header=None)
bootres_bootsamples = pd.read_csv(f"{results_dir}/bootres_bootsamples.csv", header=None)
bootres_x_weights_normed = pd.read_csv(f"{results_dir}/bootres_x_weights_normed.csv", header=None)
bootres_x_weights_stderr = pd.read_csv(f"{results_dir}/bootres_x_weights_stderr.csv", header=None)
permres_permsamples = pd.read_csv(f"{results_dir}/permres_permsamples.csv", header=None)
permres_pvals= pd.read_csv(f"{results_dir}/permres_pvals.csv", header=None)
varexp = pd.read_csv(f"{results_dir}/varexp.csv", header=None)
x_weights= pd.read_csv(f"{results_dir}/x_weights.csv", header=None)
y_loadings= pd.read_csv(f"{results_dir}/y_loadings.csv", header=None)
splitres_u_vcorr_pvals = pd.read_csv(f"{results_dir}splitres/splitres_u-vcorr_pvals.csv", header=None)


#%% Visualization: P-values and covariance explained by each latent variable

# Initialize figure
fig, ax1 = plt.subplots()
LV = list(range(1,df.shape[1]+1))

# Covariance explained in red
color = 'tab:red'
ax1.set_xlabel('LV')

ax1.set_ylabel('% covariance', color=color)
ax1.scatter(LV, varexp*100, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.tick_params(axis='x')

# Instantiate a second axes that shares the same x-axis
ax2 = ax1.twinx()  

# P-values in blue
color = 'tab:blue'
ax2.set_ylabel('p-value', color=color)
ax2.scatter(LV, permres_pvals, color=color)
ax2.tick_params(axis='y', labelcolor=color)

# Add horizontal blue line at 0.05
plt.axhline(y=0.05, color="tab:blue")
ax1.xaxis.grid(which='major')

fig.tight_layout()

# Title and font sizes
ax1.set_title("Covariance explained and p-values")

ax1.title.set_fontsize(20)
ax1.xaxis.label.set_fontsize(20)
ax1.yaxis.label.set_fontsize(20)
ax2.yaxis.label.set_fontsize(20)
ax1.tick_params(axis='x', labelsize=15)
ax1.tick_params(axis='y', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])

# Save figure
fig.savefig(f"{vis_dir}/pval_by_covexp_per_LV.png",bbox_inches = "tight", dpi=1000)

#%% Visualization: Split-half resampling p-values

# Initialize figure
fig, ax1 = plt.subplots()
LV = list(range(1,df.shape[1]+1))

LV_sp1=list()
LV_sp2=list()
LV_sp1[:] = [i-0.1 for i in LV]
LV_sp2[:] = [i+0.1 for i in LV]

# Ucorr in red
color = 'tab:red'
ax1.set_xlabel('LV')

ax1.set_ylabel('P Ucorr', color=color)
ax1.scatter(LV_sp1, np.minimum(splitres_u_vcorr_pvals[0], 0.199), color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.tick_params()
ax1.set_ylim([0, 0.20])

ax2 = ax1.twinx()

# Vcorr in blue
color = 'tab:blue'
ax2.set_ylabel('P Vcorr', color=color)
ax2.scatter(LV_sp2, np.minimum(splitres_u_vcorr_pvals[1], 0.199), color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim([0, 0.20])

# Add horizontal blue line at 0.05
plt.axhline(y=0.05, color="tab:blue")
ax1.xaxis.grid(which='major')

fig.tight_layout()

# Title and fonts
ax1.set_title("Split-half resampling p-values")

ax1.title.set_fontsize(20)
ax1.xaxis.label.set_fontsize(20)
ax1.yaxis.label.set_fontsize(20)
ax2.yaxis.label.set_fontsize(20)
ax1.tick_params(axis='x', labelsize=15)
ax1.tick_params(axis='y', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])

# Save figure
fig.savefig(f"{vis_dir}/splitres_u_vcorr_pvals.png",bbox_inches = "tight", dpi=1000)


#%% Visualization: behavioral loadings on LVs

behav = list(Y_raw.columns)
y_pos = list(range(0,len(behav)))

# Plot every significant LV
for i in range(permres_pvals.shape[0]): # i ranges from 0 to (num_LVs – 1)
    if (permres_pvals[0][i] <= 0.05):  
        y_load = y_loadings.iloc[:,i]
        x_error = []
        
        # Load data
        for j in range(len(behav)):
            y_loadings_ci = np.loadtxt(f"{results_dir}y_loadings_ci/bootres_y_loadings_ci_behaviour_{j}.csv",delimiter=",")
            pair = y_loadings_ci[i]
            entry = [math.fabs(float((y_load[j]-pair[0]))),math.fabs(float((y_load[j]-pair[1])))]
            x_error.append(entry)
            
        x_error = np.array(x_error)
        x_error = x_error.T
        
        # Assign colors for significant variable loadings
        colors=[]
        for k in range(x_error.shape[1]):
            a=y_load[k]-x_error[0,k]
            b=y_load[k]+x_error[1,k]
            c=a*b
            if (c >= 0):
                colors.append("green")
            elif (c < 0):
                colors.append("grey")
        
        # Make plot
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
        plt.savefig(f"{vis_dir}/LV_{i+1}_behav", bbox_inches='tight')
        plt.show()    

#%% Visualization: Brain loadings on LVs

# Plot every significant LV
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
    
        BSR_thresh= 3.29 #equivalent to p=0.001

        # Assign colors for significant variable loadings
        colors=[]
        for k in range(len(y_all['y_load'])):
            c=y_all['y_load'][k]
            if (abs(c) >= BSR_thresh):
                if ('CT' in y_all['brain'][k]):
                    colors.append("red")
                elif ('TBV' in y_all['brain'][k]):
                    colors.append("red")
                elif ('ICV' in y_all['brain'][k]):
                    colors.append("red")
                else:
                    colors.append("blue")
                    
            elif (abs(c) < BSR_thresh):
                colors.append("grey")
        
        # Make plot
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
        plt.axvline(x=BSR_thresh, color="black")
        plt.axvline(x=-BSR_thresh, color="black")
        ax.yaxis.grid()
        plt.savefig(f"{vis_dir}/LV_{i+1}_brain", bbox_inches='tight')
        plt.show()    