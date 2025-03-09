import pandas as pd
import pulp
import os

os.chdir('/home/linuxbrew/workspace')

print("Load sra data file")
filtered_data = pd.read_excel('filtered.xlsx', index_col=0)

antibiotic_list = ['ciprofloxacin', 'cefotaxime', 'ceftazidime', 'gentamicin'] # edit this list for selected antibiotic
minimum_resistant = 1000 # edit this value for minimum resistant genome
minimum_susceptible = 1000 # edit this value for minimum susceptible genome

sample_list = filtered_data.index.to_list()
# create dict for all antibiotic: with each antibiotic antibiotic_r for resistant, antibiotic_s for susceptible, save only sra accession as str list
antibiotic_dict = {}
for antibiotic in antibiotic_list:
  antibiotic_dict[antibiotic + '_r'] = filtered_data[filtered_data[antibiotic] == 1].index.to_list()
  antibiotic_dict[antibiotic + '_s'] = filtered_data[filtered_data[antibiotic] == 0].index.to_list()

# create LP problem select at least resistant and susceptible sra for each antibiotic and minimize the total number of selected sra
prob = pulp.LpProblem("SelectSRA", pulp.LpMinimize)
# create variables for each sra
sra_vars = pulp.LpVariable.dicts("SRA", sample_list, 0, 1, pulp.LpBinary)
# create objective function
prob += pulp.lpSum([sra_vars[sra] for sra in sample_list])
# create constraints
for antibiotic in antibiotic_list:
  prob += pulp.lpSum([sra_vars[sra] for sra in antibiotic_dict[antibiotic + '_r']]) >= min(minimum_resistant, len(antibiotic_dict[antibiotic + '_r']))
  prob += pulp.lpSum([sra_vars[sra] for sra in antibiotic_dict[antibiotic + '_s']]) >= min(minimum_susceptible, len(antibiotic_dict[antibiotic + '_s']))

prob.solve()
selected_sra = [sra for sra in sample_list if sra_vars[sra].value() == 1]

# save selected sra to file
with open('selected_sra.txt', 'w') as f:
  for sra in selected_sra:
    f.write(sra + '\n')

# select data to excel
filtered_data.loc[selected_sra,antibiotic_list].to_excel('selected_data.xlsx')


