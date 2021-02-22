from sys import argv
script, in_file, guuids, out = argv
fasta= open(in_file)
newnames= open(guuids)
newfasta= open(out, 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write('>'+newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()
