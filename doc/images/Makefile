all: tree.png nodetable.png
	dot -Kfdp -n -Tpng tree.dot -o tree.png -Gsize=3,5\!
	dot -Kfdp -n -Tpng nodetable.dot -o nodetable.png -Gsize=3,5\!

nodetable.png: nodetable.dot

