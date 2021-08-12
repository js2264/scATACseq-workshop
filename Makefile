.PHONY: run_local build

build: 
	find content/ -name *html -delete
	Rscript -e "blogdown::build_dir('content/')"
	Rscript -e "blogdown::build_dir('static/Exercises/')"
	Rscript -e "blogdown::build_dir('static/Presentations/')"
	Rscript -e "blogdown::build_site()"

run_local:
	Rscript -e "blogdown::hugo_server(host='127.0.0.1', port = 4321)"

