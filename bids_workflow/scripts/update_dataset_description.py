import json

with open(snakemake.input.json_file,'r') as f:
    json_dict = json.load(f)

json_dict['Name'] = snakemake.wildcards.prefix

with open(snakemake.output.json_file,'w') as f:
    json.dump(json_dict, f, indent=4)


