all: tree.png tables.png
	dot -Kfdp -n -Tpng tree.dot -o tree.png -Gsize=3,5\!
	dot -Kfdp -n -Tpng tables.dot -o tables.png -Gsize=3,5\!

tables.png: tables.dot

