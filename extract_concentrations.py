import os
import sys
import csv

from chemdataextractor import Document



if len(sys.argv) == 1:
	print('No file path provided. Terminating program.')
	sys.exit()


file_path = sys.argv[1]
if not os.path.isfile(file_path):
	print('Path does not refer to a valid file. Terminating program')
	sys.exit()

with open(file_path, 'rb') as f:
	#Convert to Document
	full_text_doc = Document.from_file(f)



#Extract lists of records from Documents
print("Extracting... (This may take a minute or two)")
doc_records = full_text_doc.records.serialize()
print("Extracted")
print(doc_records)

concentration_matrix = []
cleaned_doc_records = [record for record in doc_records if 'measured_concentrations' in record]
for record in cleaned_doc_records:
	name = record['names'][0]
	for concentration in record['measured_concentrations']:
		concentration_matrix.append([name, *concentration.values()])


file_name = os.path.basename(file_path).split('.')[0]#Get passed filename without suffix
with open('concentrations_{}.tsv'.format(file_name), 'w', newline='', encoding='utf-8') as tsvfile:
	writer = csv.writer(tsvfile, delimiter='\t')
	writer.writerow(['Name', 'Value', 'Units'])
	print(concentration_matrix)
	writer.writerows(concentration_matrix)