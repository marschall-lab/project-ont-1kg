import sys

for line in sys.stdin:
    if line.startswith('#'):
        print(line.strip())
        continue
    fields = line.split()
    if fields[6] == '.':
        fields[6] = 'PASS'
    print('\t'.join(fields))