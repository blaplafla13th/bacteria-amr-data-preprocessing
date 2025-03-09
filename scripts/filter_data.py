import os

import pandas as pd

os.chdir('/home/linuxbrew/workspace')

print("Load genome file")
genome = pd.read_csv('BVBRC_genome.csv', header=0, index_col=0)
genome.info()
print("Filter genome file")
# Filter columns, you can edit from here
genome_filtered = genome.loc[:,
                  ['Genome Name', 'Genome Status', 'Genome Quality', 'SRA Accession', 'Sequencing Platform']]
# Drop platform null
genome_filtered = genome_filtered.dropna(subset=['Sequencing Platform'])
# Platform is containing illumina, hiseq nextseq novaseq truseq miseq, ignore case
genome_filtered = genome_filtered[
  genome_filtered['Sequencing Platform'].str.contains('illumina|hiseq|nextseq|novaseq|truseq|miseq', case=False)]
# Genome quality Good only
genome_filtered = genome_filtered[genome_filtered['Genome Quality'] == 'Good']
# Drop SRA null
genome_filtered = genome_filtered.dropna(subset=['SRA Accession'])
# split SRA Accession by ; then take the first one have letter R at third position
genome_filtered = genome_filtered[genome_filtered['SRA Accession'].str.split(',').str[0].str[2] == 'R']
genome_filtered['SRA Accession'] = genome_filtered['SRA Accession'].str.split(',').str[0]
# End of editing
genome_filtered.info()

print("Load AMR file")
amr = pd.read_csv('BVBRC_genome_amr.csv', header=0, index_col=1)
amr.info()
print("Filter AMR file")
# Filter columns, you can edit from here
amr_filtered = amr.loc[:, ['Antibiotic', 'Resistant Phenotype']]
# AMR filter Resistance Phenotype: Susceptible and Resistant only
amr_filtered = amr_filtered[amr_filtered['Resistant Phenotype'].isin(['Susceptible', 'Resistant'])]
# End of editing
amr_filtered.info()

print("Join genome and AMR data")
data = genome_filtered.join(amr_filtered, how='inner')
data = data.reset_index()
# create new data frame have index is SRA Accession, and have columns Genome ID, Genome Status, Genome Quality, Sequencing Platform, and create new columns for each antibiotic, if the genome is resistant to the antibiotic, the value is 1, otherwise 0
antibiotics = data['Antibiotic'].unique()
antibiotic_values = {antibiotic: 0 for antibiotic in antibiotics}
data_dict = {}
for row in data.iterrows():
  index = row[1]['SRA Accession']
  if index not in data_dict:
    data_dict[index] = {'Genome ID': row[1]['Genome ID'], 'Genome Status': row[1]['Genome Status'],
                        'Genome Quality': row[1]['Genome Quality'],
                        'Sequencing Platform': row[1]['Sequencing Platform']}
    for antibiotic in antibiotics:
      data_dict[index][antibiotic] = -1
  data_dict[index][row[1]['Antibiotic']] = 1 if (
      row[1]['Resistant Phenotype'] == 'Resistant') else 0 if row[1]['Resistant Phenotype'] == 'Susceptible' else -1
  if data_dict[index][row[1]['Antibiotic']] == 1:
    antibiotic_values[row[1]['Antibiotic']] += 1
# sort antibiotic_values by value
antibiotic_values = dict(sorted(antibiotic_values.items(), key=lambda item: item[1], reverse=True))
# create data frame from data_dict column order Genome ID, Genome Status, Genome Quality, Sequencing Platform,
# and antibiotics in sorted order in antibiotic_values
data_new = pd.DataFrame(data_dict).T
data_new = data_new[
  ['Genome ID', 'Genome Status', 'Genome Quality', 'Sequencing Platform'] + [antibiotic for antibiotic in
                                                                             antibiotic_values]]
# sort index
data_new = data_new.sort_index()
data_new.info()

print("Save filtered data")
data_new.to_excel('filtered.xlsx')
