---
title: "3. Rstudio"
---

Rstudio offers a graphical interface 
to facilitate the interaction between a user and an underlying 
programming language (this is sometimes called IDE, or 
*integrated development environment*). It can be very useful when a user is 
not necessarily proficent with command line-based computing. 
However, such graphical interfaces are not 
always able to connect to services such as AWS.  

Since most of the analyses we will conduct will be done on AWS, 
we'd like to be able to use Rtudio directly from there. 
To access RStudio running on the AWS instance from the web, 
simply go to the following address: 

```sh
"http://${IP}:8787"
```

Don't forget, you can run a `bash` terminal from within RStudio! 
This may come handy if you want to process some data with `cellranger`, for example. 
To do this, simply click on the `terminal` button next on the bottom left panel. 
You should now be in your own `${home}` directory (`~`). 
