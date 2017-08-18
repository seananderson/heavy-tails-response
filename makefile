all: anderson-et-response.pdf anderson-et-response.docx

anderson-et-response.pdf: anderson-et-response.md refs.bib pnas.csl
	pandoc anderson-et-response.md --bibliography=refs.bib --csl=pnas.csl -o anderson-et-response.pdf

anderson-et-response.docx: anderson-et-response.md refs.bib pnas.csl
	pandoc anderson-et-response.md --bibliography=refs.bib --csl=pnas.csl -o anderson-et-response.docx
