import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
plt.rc('xtick',labelsize=6)
plt.rc('ytick',labelsize=8)

datasets = pd.read_table("./data/dataset_config_test.csv")
yeast_ref = pd.read_table("./data/yeast_reference.csv")


sigs = []
groups = []
sims = []
for i in range(datasets.shape[0]):
    sigs.append(pd.read_csv("./data/"+datasets.iloc[i,1],low_memory=False,sep='\t'))
    groups.append(pd.read_csv("./data/"+datasets.iloc[i,2],sep='\t', dtype={"start":int,"end":int,"group_id":int}))
    sims.append(pd.read_csv("./data/"+datasets.iloc[i,3],sep='\t', header=None))
    
for i in range(datasets.shape[0]):
    sigs[i].iloc[:,3:] = sigs[i].iloc[:,3:].fillna(0)
    
for i in range(datasets.shape[0]):
    sigs[i].iloc[:,1:3] = sigs[i].iloc[:,1:3].fillna("N/A")


#Get group 
def getgroup(dataset):
    return groups[dataset]

#Find index by idr name
def getindex_idr(name, dataset):
    test = sigs[dataset].iloc[:,0]==name
    if test.any()  == False:
        return -1
    else:
        return sigs[dataset].index.values[sigs[dataset].iloc[:,0]==name][0]  

# Find idrname by protein common name
def getindex_cid(name, dataset):
    index = []
    for i in range(len(sigs[dataset])):
        cid =  sigs[dataset].iloc[i,2]  # in this table, third column is common name separated by ";"
        if  name.lower() == cid.lower():
            index.append(sigs[dataset].iloc[i,0])
    return index

#Find idrname by systematic name
def getindex_sid(name, dataset):
    index = []
    for i in range(len(sigs[dataset])):
        sid =  sigs[dataset].iloc[i,1]
        if sid.lower() == name.lower():
            index.append(sigs[dataset].iloc[i,0])
    return index

# Find the list of similar IDRs, the first item is the orginal IDR
def getsimi(IDR, dataset):
    return sims[dataset].iloc[IDR,:].tolist()

# Find common name from index
def getname(index, dataset):
    return sigs[dataset].iloc[index,2]

# Find sysmatic name/uniprot id from index
def getsys(index, dataset):
    return sigs[dataset].iloc[index,1]

# Find Yeast uniprot id from sysmatic name
def getuni(sys):
    for i in range(len(yeast_ref)):
        if yeast_ref.iloc[i,0] == sys:
            return yeast_ref.iloc[i,1]
    return None   

def get_dataset_name(index):
    return datasets.iloc[index,0]

# Data figures
def sigviz(IDR, Format="bar", Group=1, dataset=0):    
    i=IDR
    
    # ZScore Values of every group sorted by absolute value
    group_sort = [None for _ in range(len(groups[dataset]))] 
    #for x in range(datasets.shape[0]):
    for y in range(len(groups[dataset])):
        group_sort[y] = sorted(sigs[dataset].iloc[i,groups[dataset].iloc[y,3]-1:groups[dataset].iloc[y,4]], key =abs, reverse=True)
        
    # Sortd absolute ZScore values for features in every group
    group_abs = [None for _ in range(len(groups[dataset]))] 
    #for x in range(datasets.shape[0]):
    for y in range(len(groups[dataset])):
        group_abs[y] = sorted(pd.Series.tolist(abs(sigs[dataset].iloc[i,groups[dataset].iloc[y,3]-1:groups[dataset].iloc[y,4]])), reverse=True)
    
    #Calculate number of features in every group
    group_len = [groups[dataset].iloc[y,4] - groups[dataset].iloc[y,3] + 1 for y in range(len(groups[dataset]))] 
     
    #Group with maximum 5 feature length
    group_len_limited = [np.minimum(group_len[y] ,5) for y in range(len(groups[dataset]))] 


    
    if Format=="bar": #Draw bar plot
        barWidth = 1
      
        fig1 = plt.figure()
        ax1 = SubplotHost(fig1, 111)
        fig1.add_subplot(ax1) 
        
        # Make the plot
        colors = ['#36072d', '#910f3f', '#9c0615', '#d45a26', '#edc92b']
        for x in range(groups[dataset].shape[0]):
            plt.bar(sum(len for len in group_len_limited[0:x]) + np.arange(np.minimum(group_len[x],5))+1, group_abs[x][0:5],
                    color=colors[x%5], width=barWidth, edgecolor='white')


     
        # Add xticks on the middle of the group bars
        plt.ylabel('z-score', fontweight='bold')      
        ax1.xaxis.set_major_locator(MultipleLocator(sum(group_len_limited)/(len(group_len_limited)-1)))
        labels = [" "] +[group for group in groups[dataset].iloc[:,1] ] + [" "]
        ax1.set_xticklabels(labels,ha='right')
 
        ax2 = ax1.twiny()
        offset = 0, -20 # Position of the second axis
        new_axisline = ax2.get_grid_helper().new_fixed_axis
        ax2.axis["bottom"] = new_axisline(loc="bottom", axes=ax2, offset=offset)
        ax2.axis["top"].set_visible(False)
        #calculate where to cut group 
        count_mean=groups[dataset].iloc[:,0]=="Mean"
        group_cut=count_mean.value_counts()[1]/groups[dataset].shape[0]
        ax2.set_xticks([0.0, group_cut, 1.0])
        ax2.xaxis.set_major_formatter(ticker.NullFormatter())
        ax2.xaxis.set_minor_locator(ticker.FixedLocator([group_cut/2, (1 - group_cut)/2 + group_cut]))
        ax2.xaxis.set_minor_formatter(ticker.FixedFormatter(['mean', 'log variance']))
    
        ax3 = ax1.twiny()
        offset = 0, -35
        new_axisline = ax3.get_grid_helper().new_fixed_axis
        ax3.axis["bottom"] = new_axisline(loc="bottom", axes=ax3, offset=offset)
        ax3.axis["top"].set_visible(False)

        ax3.set_xticks([0.0, 1.0])
        ax3.xaxis.set_major_formatter(ticker.NullFormatter())
        ax3.xaxis.set_minor_locator(ticker.FixedLocator([0.5]))
        ax3.xaxis.set_minor_formatter(ticker.FixedFormatter(['Molecular Features']))
    
        
        figbar = plt.gcf()
        figbar.subplots_adjust(bottom=0.22)
        figbar.set_size_inches(6, 3)
        figbar.savefig('./static/image/sigbar'+str(i)+'.png', dpi=100)
 
        return()
        
    if Format == "div": #Draw diversion plot
        mean_or_variance = groups[dataset].iloc[Group-1,0]
        group_name = groups[dataset].iloc[Group-1,2]
        colnames = sigs[dataset].columns.values
        feature = colnames[groups[dataset].iloc[Group-1,3]-1:groups[dataset].iloc[Group-1,4]]
        score_list = group_sort[Group-1][0:10]
        score_list=[x for x in reversed(score_list) if ~np.isnan(x)]
        fvalue=sigs[dataset][feature].iloc[i]
        feature_list=[None]*len(score_list)
        for x1 in range(len(score_list)):
            for j in range(len(fvalue)):
              if (score_list[x1]==fvalue[j]):
                  feature_list[x1]=fvalue.index[j]
                  
        fig, yx = plt.subplots(constrained_layout=True)
        plt.title('Top {0} Z-scores in \n {1} Group '.format(mean_or_variance, group_name), 
                  fontdict={'size':10})
        plt.hlines(np.arange(len(score_list)), xmin=0, xmax=score_list)
        for x, y, tex in zip(score_list, np.arange(len(score_list)), score_list):
            plt.text(x, y, round(tex, 2), horizontalalignment='right' if x < 0 else 'left', 
                         verticalalignment='center', fontdict={'color':'red' if x < 0 else 'green', 'size':6})
        plt.yticks(np.arange(len(score_list)), feature_list, fontsize=6)
        plt.grid(linestyle='--', alpha=0.35)
        plt.gca().set(ylabel='$Molecular-features$', xlabel='$Z-scores$')
        
        figdiv = plt.gcf()
        figdiv.set_size_inches(5, 3.5)
        figdiv.savefig('./static/image/sigdiv'+str(i)+'g'+str(Group)+'.png', dpi=100)
        plt.clf()
        
        return()


#Draw profile plot
def sigpro(IDR, dataset):

    values = sigs[dataset].iloc[IDR,3:]
    fig1 = plt.figure()
    fig1.set_size_inches(7, 3.5)
    plt.subplots_adjust(bottom=0.20)
    ax1 = SubplotHost(fig1, 111)
    fig1.add_subplot(ax1) 
    plt.ylabel("Z-Scores")
    ax1.stem(values, markerfmt=' ',use_line_collection=True,basefmt='k-')
    ax1.xaxis.set_major_locator(MultipleLocator(len(values)/(groups[dataset].shape[0]-1)))
    labels = [" "] +[group for group in groups[dataset].iloc[:,1] ]
    ax1.set_xticklabels(labels, horizontalalignment="right")
    plt.grid(axis='y')

    max=np.amax(values)
    max_x=np.argmax(np.array(values))
    min=np.amin(values)
    min_x=np.argmin(np.array(values))
    ax1.annotate(sigs[dataset].columns[max_x+3],fontsize=7,
            xy=(max_x, max), xycoords='data',bbox=dict(boxstyle="round", fc="0.8"),
            xytext=(15, 15), textcoords='offset points',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="angle,angleA=0,angleB=90,rad=10"))

    ax1.annotate(sigs[dataset].columns[min_x+3],fontsize=7,xy=(min_x, min), xycoords='data',
            xytext=(12, 12), textcoords='offset points',
            bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="angle,angleA=90,angleB= 0,rad=10"))
    for k in range(4, int(plt.ylim()[1])):
        ax1.axhspan(k, k+1, color='sandybrown', alpha=0.3)
    for k in range(int(plt.ylim()[0]),-4):
        ax1.axhspan(k, k+1, color='sandybrown', alpha=0.3)
 
    ax2 = ax1.twiny()
    offset = 0, -20 # Position of the second axis
    new_axisline = ax2.get_grid_helper().new_fixed_axis
    ax2.axis["bottom"] = new_axisline(loc="bottom", axes=ax2, offset=offset)
    ax2.axis["top"].set_visible(False)
    #calculate where to cut group 
    count_mean=groups[dataset].iloc[:,0]=="Mean"
    group_cut=count_mean.value_counts()[1]/groups[dataset].shape[0]
    ax2.set_xticks([0.0, group_cut, 1.0])
    ax2.xaxis.set_major_formatter(ticker.NullFormatter())
    ax2.xaxis.set_minor_locator(ticker.FixedLocator([group_cut/2, (1 - group_cut)/2 + group_cut]))
    ax2.xaxis.set_minor_formatter(ticker.FixedFormatter(['mean', 'log variance']))
    
    ax3 = ax1.twiny()
    offset = 0, -35
    new_axisline = ax3.get_grid_helper().new_fixed_axis
    ax3.axis["bottom"] = new_axisline(loc="bottom", axes=ax3, offset=offset)
    ax3.axis["top"].set_visible(False)

    ax3.set_xticks([0.0, 1.0])
    ax3.xaxis.set_major_formatter(ticker.NullFormatter())
    ax3.xaxis.set_minor_locator(ticker.FixedLocator([0.5]))
    ax3.xaxis.set_minor_formatter(ticker.FixedFormatter(['Molecular Features']))
    
    figpro = plt.gcf()
    figpro.set_size_inches(6, 3)
    figpro.savefig('./static/image/sigpro'+str(IDR)+'.png', dpi=100)
    plt.clf()










