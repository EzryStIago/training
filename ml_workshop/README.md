# Introduction

This is an introductory overview of machine learning. We will not dive into code, though some code samples are available. One valuable part of the workshop is that participants share what kind of data they are working on. 

## Here are the resources we are going to be using or referencing: 
- slides: https://docs.google.com/presentation/d/1b12AuZZlJRUx_22ojqImKIQ1WPaWtmG5oNvV7IZHAyE/edit?usp=sharing 
- https://datastudio.google.com         - create dashboards and visualizations, like Tableau
- https://colab.research.google.com  - free GPU in Jupyter-like environment
- https://ondemand.hpc.rutgers.edu  - Jupyter notebook on Rutgers’ Amarel cluster
- https://jupyterhub.hpc.rutgers.edu  - additional place where you can pull workshop data (VPN)

## Jupyter shortcuts

- 2 types of cells - code and text
- 2 modes of interacting with cells - “editing” and “command”
- CTRL-Enter for executing the cell
- ESC to enter the command mode
- m to change the cell to text (i.e. “markdown”)
- CTRL-/  for commenting out a piece of code
- dd for deleting the cell
- a for new cell above, b for new cell below
- CTRL-minus for splitting a cell

## Python most important tips: 

- tab-completion for available functions on the object
- ? gives the documentation on that object or function


## Datasets used in the workshop (in data folder): 

- Titanic - https://www.kaggle.com/c/titanic 
- Breast cancer - included with scikit-learn
- Adult income - https://archive.ics.uci.edu/ml/datasets/adult 

## On campus help:

- OARC Data Science office hours - Mondays 3-4:30pm in CoRE 710 - please email ahead
- Graduate student specialist program - Library


## Additional selected tutorials and courses: 

- https://www.fast.ai/  by Jeremy Howard - excellent practical courses - python
- https://www.youtube.com/playlist?list=PLOg0ngHtcqbPTlZzRHA2ocQZqB1D_qZ5V “Introduction to Statistical Learning Series” by Hastie and Tibshirani (Stanford) - R - has an accompanying book “An Introduction to Statistical Learning: with Applications in R” 

# Technical instructions how to access materials: 

- **slides**: https://docs.google.com/presentation/d/1b12AuZZlJRUx_22ojqImKIQ1WPaWtmG5oNvV7IZHAyE/edit?usp=sharing
- **jupyter notebook** (available options): 
  * download from: https://jupyterhub.hpc.rutgers.edu/kp807/titanic_notebook_clean_master.ipynb 
  * copy on Amarel from /projects/oarc/users/share/kp807 to your home directory
  * save to your google drive from this shared link: https://colab.research.google.com/drive/12V3yfjPfkHvE8SfCqINVgBd_S07PGImS
- **data** (available options): 
  * download from: https://jupyterhub.hpc.rutgers.edu/kp807/   (train.csv and test.csv) and upload
  * copy on Amarel from /projects/oarc/users/share/kp807 to your home directory
  * save to your google drive from this shared link: https://datastudio.google.com/open/1m1T50t6qbaw7y7xwuVQLEo-K5zgfEBEo

# Google Cloud options 

- you must have a google account (such as your rutgers scarletmail account) - every gmail account gets $300 free cloud credits for 1 year
- exploring data and creating dashboards: https://datastudio.google.com/ (max size of 100Mb for exploration)
- machine learning on tabular datasets: https://cloud.google.com/automl-tables/  (click on "Go to console")
- Google Colab (hosted Jupyter notebook): https://colab.research.google.com/

# Amarel instructions

- https://ondemand.hpc.rutgers.edu -> Interactive Apps -> Comp. Gen. Jupyter  (with partition main- TODO:change)
- copying data and notebook to your directory: 
```
mkdir mlworkshop
cp /projects/oarc/users/share/kp807/*csv /home/$USER/mlworkshop
cp /projects/oarc/users/share/kp807/*ipynb /home/$USER/mlworkshop
```

