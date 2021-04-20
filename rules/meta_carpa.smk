# Meta-analysis using METACARPA https://github.com/hmgu-itg/metacarpa/

rule metacarpa_download:
	input: HTTP.remote("https://github.com/hmgu-itg/metacarpa/releases/download/1.0.1/metacarpa", keep_local=True)
	output: "resources/metacarpa/metacarpa"