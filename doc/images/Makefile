all: tree.png tables.png structures.png
	dot -Kfdp -n -Tpng tree.dot -o tree.png -Gsize=3,5\!
	dot -Kfdp -n -Tpng tables.dot -o tables.png -Gsize=3,5\!
	dot -n -Tpng structures.dot -o structures.png 

tables.png: tables.dot
structures.png: structures.dot

