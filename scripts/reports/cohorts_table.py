import pandas as pd

# open case/control count spreadsheet
casecontrol= pd.read_excel(snakemake.input[0])

# clean up cohorts filename text
casecontrol['Dataset'] = casecontrol['Dataset'].map(lambda dataset: dataset.replace('.aligned.gz', ""))

# write out as a plain table
with open(snakemake.output[0], 'w') as out:
    casecontrol.to_string(out)