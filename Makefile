update-readme:
	export PATH="/Users/davide/Library/Python/3.9/bin:$$PATH"
	jupyter nbconvert --to markdown lic_examples.ipynb
	mv lic_examples.md README.md
