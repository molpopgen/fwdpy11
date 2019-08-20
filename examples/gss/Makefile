all: pop.lzma stats.sqlite3 stats.pdf compare_genetic_values.png mean_trait_value_from_tree_sequences.png

clean:
	rm -f pop.lzma *.sqlite3 *.pdf *.png 

pop.lzma: DiploidPopulationGSSmo.py 
	python3 DiploidPopulationGSSmo.py --popsize 1000 --mu 5e-3 --sigma 0.25 --filename pop.lzma --preserve 1 --num_ind 1000 --record --statfile stats.sqlite3 

stats.sqlite3: pop.lzma DiploidPopulationGSSmo.py

stats.pdf: pop.lzma plotstats.R
	./plotstats.R stats.sqlite3 9950 10101

mean_trait_value_from_tree_sequences.png: pop.lzma plot_genetic_values_from_tree_sequences.py
	python3 plot_genetic_values_from_tree_sequences.py pop.lzma

compare_genetic_values.png: pop.lzma iterate_variants_in_tree_sequences.py
	python3 iterate_variants_in_tree_sequences.py pop.lzma
