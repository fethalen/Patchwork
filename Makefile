format:
	jlfmt --threads=auto -- --inplace .

format-check:
	jlfmt --check -v .
