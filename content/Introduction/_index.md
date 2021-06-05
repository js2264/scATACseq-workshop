---
title: " Single-cell RNA-seq analysis with R/Bioconductor"
---

## Course content: 

The workshop content is available [at this adress: https://github.com/js2264/...](https://github.com/js2264/...). 

You can clone it locally with `git`: 

```sh
git clone https://github.com/js2264/....git
```

## To log in RStudio server: 

[Simply go to the following web page](http://54.188.62.211:8787/): `http://54.188.62.211:8787/`

## IP for today (2021-01-14): 

```sh
IP="54.188.62.211"
```

## To log in your AWS session: 

```sh
KEY="path/to/key.pem"
chmod 600 "${KEY}"
USER="myuserID"
IP="54.188.62.211"
ssh -i "${KEY}" "${USER}"@"${IP}"
```

More information is available [here](/{{<myPackageUrl>}}how_to_log_in_aws/#connect-to-pre-configured-aws-machine) if needed. 
