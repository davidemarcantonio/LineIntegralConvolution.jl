update-readme:
	( \
		source env/bin/activate; \
		jupyter nbconvert --to markdown lic_examples.ipynb; \
		mv lic_examples.md README.md; \
		deactivate; \
	)

julia:
	sudo julia

# configure-venv:
# 	python3 -m venv env
# 	source env/bin/activate
# 	pip install --upgrade pip
# 	pip install jupyterlab
# 	deactivate

# delete-venv:
# 	rm -rf env
	