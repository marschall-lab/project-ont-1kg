import sys
from collections import defaultdict

# list of IDs to skip
var_ids_file = sys.argv[1]
var_ids = defaultdict(lambda: False)

for line in open(var_ids_file, 'r'):
	var_ids[line.strip()] = True

sys.stderr.write('Read ' + str(len(var_ids)) + ' IDs from input list.\n')

total_lines = 0
skipped_lines = 0
for line in sys.stdin:
	if line.startswith('#'):
		print(line[:-1])
		continue
	fields = line.split()
	assert fields[2] != '.'
	total_lines += 1
	if var_ids[fields[2]]:
		skipped_lines += 1
		continue
	print(line[:-1])

sys.stderr.write('Skipped ' + str(skipped_lines) + ' of ' + str(total_lines) + ' from input vcf.\n')	
